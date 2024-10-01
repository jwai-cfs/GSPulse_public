import tkinter as tk
from tkinter import ttk
from functools import partial
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
import matplotlib.pyplot as plt
from shape_callbacks import shape_create_deadstart, interparc
from intersections import intersection
import numpy as np
import json
import argparse
from tok_geo_defaults import tok_geo_defaults
from pathlib import Path
import os

class ShapeApp:
    """
    CLASS: ShapeApp
    DESCRIPTION: 
    """

    def __init__(self, tokamak):
        """
        METHOD: __init__
        DESCRIPTION:                
        """        
        
        # define root window
        self.root = tk.Tk()
        self.define_root_window()                              

        # initialize tokamak and defaults
        self.tokamak = tokamak        
        self.geo = tok_geo_defaults(tokamak)                
        shp = {}
        shp = self.add_missing_shape_fields(shp)
        shp = self.add_aux_geom_params(shp)
        self.tk_shape = self.dict2tkdict(shp)
        self.default_shape_fn = f'{str(Path(__file__).parent)}/default_shapes/{self.tokamak}_shape.json'        
        if os.path.exists(self.default_shape_fn):
            self.load_shape_default()           
                             
        # create notebook with tabs
        notebook = ttk.Notebook(self.root)   
        tab1 = ttk.Frame(notebook)        
        tab2 = ttk.Frame(notebook)
        notebook.add(tab1, text='shape')
        notebook.add(tab2, text='weights')
        notebook.pack(expand=1, fill='both')

        # fileio panel
        self.add_fileio_panel(tab1)
        
        # add panels for shape parameter inputs        
        self.add_shape_params_panel(tab1)
        self.add_shape_points_panel(tab1)
        self.add_segs_panel(tab1)
        self.add_plot_opts_panel(tab1)

        # add plot axes
        plot_frame = tk.Frame(tab1)
        plot_frame.pack(side='left', anchor='nw', padx=10)
        self.add_plot_axes(plot_frame)    

        # plot shape
        self.update_plots()    


    def add_missing_shape_fields(self, shape):
               
        self.n_manual_control_pts = 8
        self.n_xpts = 4
        self.n_manual_segs = 4
        
        keys =  ['rx' + str(i+1) for i in range(self.n_xpts)]
        keys += ['r' + str(i+1) for i in range(self.n_manual_control_pts)] 
        keys += ['zx' + str(i+1) for i in range(self.n_xpts)]
        keys += ['z' + str(i+1) for i in range(self.n_manual_control_pts)] 
        keys += ['Zup', 'Zlo','Rout', 'Rin','triu','tril','squo','squi','sqlo','sqli','c_xplo','c_xpup']                
        keys += ['nsegs', 'seglength', 'theta0', 'rc', 'zc', 'a', 'b']        
        for i in range(self.n_manual_segs):
            keys += [f'seg{i}_R0', f'seg{i}_Z0', f'seg{i}_Rf', f'seg{i}_Zf']

        for key in keys:
            if key not in shape.keys():
                shape[key] = np.nan
                             
        return shape       
       
    def add_aux_geom_params(self, shape):

        # auxiliary geometric parameters             
        shape['aminor']  = (shape['Rout'] - shape['Rin']) / 2.0
        shape['R0'] = (shape['Rout'] + shape['Rin']) / 2.0
        shape['bminor']  = (shape['Zup'] - shape['Zlo']) / 2.0
        shape['Z0'] = (shape['Zup'] + shape['Zlo']) / 2.0
        shape['k'] = shape['bminor'] / shape['aminor']       

        return shape

    # panel to hold fileio buttons
    def add_fileio_panel(self, parent):

        panel = tk.LabelFrame(parent, bd=0)
        panel.pack(side='bottom', anchor='sw', padx=10, pady=10)
              
        B = tk.Button(panel, text='Save Shape', command=self.save_shape)
        B.pack(side='left', anchor='sw', padx=10, pady=10)

        B = tk.Button(panel, text='Load Shape', command=self.load_shape_and_plot)
        B.pack(side='left', anchor='sw', padx=10, pady=10)
        
        B = tk.Button(panel, text='Load default', command=self.load_shape_default_and_plot)
        B.pack(side='left', anchor='sw', padx=10, pady=10)
        
        B = tk.Button(panel, text='Set as default', command=self.save_shape_default)
        B.pack(side='left', anchor='sw', padx=10, pady=10)
        
   
        
    def add_plot_opts_panel(self, parent):

        # panel to hold plot options
        panel = tk.LabelFrame(parent, text='Plot options', highlightbackground="gray", highlightthickness=2)
        panel.pack(side='bottom', anchor='nw', padx=10, pady=10)

        self.plot_segs = tk.IntVar(value=True)                 
        B = tk.Checkbutton(panel, text='plot segs', variable = self.plot_segs, command=self.update_plots)
        B.pack(side='top', anchor='nw', padx=10, pady=0) 
        
        self.label_control_pts = tk.IntVar()                 
        B = tk.Checkbutton(panel, text='label control points', variable = self.label_control_pts, command=self.update_plots)
        B.pack(side='top', anchor='nw', padx=10, pady=0)
        
        self.label_manual_control_pts = tk.IntVar()                 
        B = tk.Checkbutton(panel, text='label manual control points', variable = self.label_manual_control_pts, command=self.update_plots)
        B.pack(side='top', anchor='nw', padx=10, pady=0)

        self.label_xpts = tk.IntVar()                 
        B = tk.Checkbutton(panel, text='label x-points', variable = self.label_xpts, command=self.update_plots)
        B.pack(side='top', anchor='nw', padx=10, pady=0)

        self.hold_shape = tk.IntVar()
        B = tk.Checkbutton(panel, text='hold shapes', variable=self.hold_shape, command=self.update_plots)
        B.pack(side='top', anchor='nw', padx=10, pady=0)               
         

    def add_segs_panel(self, parent):
         
        # panel to hold segment parameter widgets
        panel = tk.LabelFrame(parent, text='Control segments', highlightbackground="gray", highlightthickness=2)
        panel.pack(side='left', anchor='nw', padx=10, pady=10)
        
        # initialize segment parameters
        seg_keys = ['nsegs', 'seglength', 'theta0', 'rc', 'zc', 'a', 'b']
        seg_key_labels = ['# segs', 'seg_length', 'theta0', 'ellipse_r0', 'ellipse_z0', 'a', 'b']
        
        # assign widgets for each segment parameter
        for i, (key, key_label) in enumerate(zip(seg_keys, seg_key_labels)):
            label = tk.Label(panel, text=key_label)
            label.grid(row=i, column=1)
            entry = tk.Entry(panel, bd=5, width=10, textvariable=self.tk_shape[key])                
            entry.bind('<Return>', self.update_plots)                                    
            entry.grid(row=i, column=3)
        
        # manual control segments
        rowstart = i
        label = tk.Label(panel, text='\n\nManual control segments')
        label.grid(row=rowstart+1, column=1, columnspan=2)
        for (i,text) in enumerate(['R0', 'Z0', 'Rf', 'Zf']):
            label = tk.Label(panel, text=text)
            label.grid(row=rowstart+2, column=i)

        rowstart += 3        
        for i in range(self.n_manual_segs):
            keys = [f'seg{i}_R0', f'seg{i}_Z0', f'seg{i}_Rf', f'seg{i}_Zf']

            for (col, key) in enumerate(keys):                
                entry = tk.Entry(panel, bd=5, width=3, textvariable=self.tk_shape[key])        
                entry.bind('<Return>', self.update_plots)                                    
                entry.grid(row=rowstart+i, column=col)

        
    def dict2tkdict(self, d):
        tkd = {}
        for key in d.keys():
            tkd[key] = tk.StringVar(value=d[key])
        return tkd

    def tkdict2dict(self, tkdict):
        d = {}
        for key in tkdict.keys():
            try:
                dum = tkdict[key].get()
                d[key] = float(dum)   # convert to numeric if possible
            except:
                d[key] = np.nan       # otherwise, nan
        return d
    
    def get_segs(self):

        # manually-defined segs
        n = self.n_manual_segs
        mansegs = np.empty((n,4))
        for i in range(n):
            keys = [f'seg{i}_R0', f'seg{i}_Z0', f'seg{i}_Rf', f'seg{i}_Zf']
            for j, key in enumerate(keys):
                dum = self.tk_shape[key].get()
                try:
                    mansegs[i,j] = float(dum)   # convert to numeric if possible
                except:
                    mansegs[i,j] = np.nan       # otherwise, nan


        # parameterized segs
        p = self.tkdict2dict(self.tk_shape)
        th = np.linspace(0, 2*np.pi, 200) + p['theta0']

        rin = p['rc'] + p['a'] * np.cos(th)
        zin = p['zc'] + p['b'] * np.sin(th)
        rout = p['rc'] + p['seglength'] * p['a'] * np.cos(th)
        zout = p['zc'] + p['seglength'] * p['b'] * np.sin(th)
        
        rout, zout = interparc(rout, zout, int(p['nsegs']))
        idx = []
        for (ro,zo) in zip(rout,zout):
            dist2 = (rin-ro)**2 + (zin-zo)**2
            idx.append(np.argmin(dist2))
        
        idx = np.asarray(idx)
        segs = np.vstack((rin[idx], zin[idx], rout, zout)).T

        # concatenate the manual and parametrized segs
        segs = np.vstack((segs, mansegs))
        return segs                       

    def plot_limiter(self, ax):
        ax.plot(self.geo['rl'], self.geo['zl'], linewidth=1.5, color='black')     

    def add_plot_axes(self, parent): 
        """
        METHOD: plot
        DESCRIPTION: 
        """       

        self.fig = Figure(figsize = (6,6), dpi = 100) 
        self.axs = [None]*3
        self.axs[0] = self.fig.add_subplot(2,2,(1,3)) 
        self.axs[1] = self.fig.add_subplot(2,2,2) 
        self.axs[2] = self.fig.add_subplot(2,2,4) 

        for i in range(3):            
            self.axs[i].grid(visible=True)            
            if i == 0:
                self.axs[i].set_ylabel('Z [m]', fontsize=12)
            if i in (0,2):
                self.axs[i].set_xlabel('R [m]', fontsize=12)
            if i == 1:
                self.axs[i].xaxis.set_ticklabels([])

        self.axs[0].set_xlim(self.geo['xlim0'])
        self.axs[0].set_ylim(self.geo['ylim0'])
        self.axs[1].set_xlim(self.geo['xlim1'])
        self.axs[1].set_ylim(self.geo['ylim1'])
        self.axs[2].set_xlim(self.geo['xlim2'])
        self.axs[2].set_ylim(self.geo['ylim2'])

        for ax in self.axs:
            ax.set_aspect('equal')

        self.fig.tight_layout()
        
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)   
        self.canvas.draw() 
        # placing the canvas on the Tkinter window 
        self.canvas.get_tk_widget().pack() 
    
        # creating the Matplotlib toolbar 
        toolbar = NavigationToolbar2Tk(self.canvas, parent) 
        toolbar.update() 
        self.canvas.get_tk_widget().pack()


   
    def add_shape_params_panel(self, parent):
        """
        METHOD: add_shape_params_panel
        DESCRIPTION:                
        """        

        # panel to hold shape parameter widgets
        shp_frame = tk.LabelFrame(parent, text='Shape parameters', highlightbackground="gray", highlightthickness=2)
        shp_frame.pack(side='left', anchor='nw', padx=10, pady=10)
            
        # initialize shape parameters
        shape_keys = ['Zup', 'Zlo','Rout', 'Rin','triu','tril','squo','squi','sqlo','sqli','c_xplo','c_xpup']
        shape_key_labels = ['Zup', 'Zlo', 'Rout', 'Rin', 'Triangularity upper', 'Triangularity lower',
                           'Squareness up/out', 'Squareness up/in', 'Squareness lo/out', 'Squareness lo/in', 'Xpt_coeff lower', 
                           'Xpt_coeff upper']
                
        # assign widgets for each shape parameter
        for i, (key, key_label) in enumerate(zip(shape_keys, shape_key_labels)):
            label = tk.Label(shp_frame, text=key_label)
            label.grid(row=i, column=2)

            entry = tk.Entry(shp_frame, 
                                bd=5, 
                                width=10, 
                                textvariable=self.tk_shape[key])                
            entry.bind('<Return>', self.update_plots)                                    
            entry.grid(row=i, column=3)

    def add_shape_points_panel(self, parent):
        """
        METHOD: add_shape_params_panel
        DESCRIPTION:                
        """        

        # panel to hold shape parameter widgets
        shp_frame = tk.LabelFrame(parent, text='Manual points', highlightbackground="gray", highlightthickness=2)
        shp_frame.pack(side='left', anchor='nw', padx=10, pady=10)

        rkeys =  ['rx' + str(i+1) for i in range(self.n_xpts)]
        rkeys += ['r' + str(i+1) for i in range(self.n_manual_control_pts)] 
        zkeys =  ['zx' + str(i+1) for i in range(self.n_xpts)]
        zkeys += ['z' + str(i+1) for i in range(self.n_manual_control_pts)] 

        # create widgets for each manual x-point and control point
        for i, (rkey, zkey) in enumerate(zip(rkeys, zkeys)):

            label = tk.Label(shp_frame, text=rkey)
            label.grid(row=i, column=1)

            entry = tk.Entry(shp_frame, 
                                bd=5, 
                                width=4, 
                                textvariable=self.tk_shape[rkey])                
            entry.bind('<Return>', self.update_plots)                                    
            entry.grid(row=i, column=2)

            label = tk.Label(shp_frame, text=zkey)
            label.grid(row=i, column=3)

            entry = tk.Entry(shp_frame, 
                                bd=5, 
                                width=4, 
                                textvariable=self.tk_shape[zkey])                
            entry.bind('<Return>', self.update_plots)                                    
            entry.grid(row=i, column=4)


    def seg_intersections(self, segs, rb, zb):
        """
        METHOD: seg_intersections
        DESCRIPTION: find intersection of control segments and boundary                
        """  
        nsegs = segs.shape[0]
        rcp = np.empty(nsegs)*np.nan
        zcp = np.empty(nsegs)*np.nan
        for i in range(segs.shape[0]):            
            rseg = segs[i,[0,2]]
            zseg = segs[i,[1,3]]
            r_, z_ = intersection(rseg, zseg, rb, zb)
            if r_.size != 0:
                rcp[i] = r_[0]
                zcp[i] = z_[0]

        return rcp, zcp                

    def set_entry_text(self, entry, text):
        """
        METHOD: set_entry_text
        DESCRIPTION:                
        """                
        entry.delete('0', 'end')
        entry.insert('0', text)


    def define_root_window(self):
        """
        METHOD: define_root_window
        DESCRIPTION:                
        """                
        self.root.title('Shape Editor')              
        
        # center the window when opened
        w = 1400                            # width for the root window
        h = 800                             # height for the root window
        ws = self.root.winfo_screenwidth()  # width of the screen
        hs = self.root.winfo_screenheight() # height of the screen   
        x = (ws/2) - (w/2)                  # calculate x and y coordinates for the window
        y = (hs/2.4) - (h/2)    
        self.root.geometry('%dx%d+%d+%d' % (w, h, x, y))              


    def update_data(self, event=None):
        
        s = self.tkdict2dict(self.tk_shape)
        s = self.add_aux_geom_params(s)        
        s['rb'], s['zb'] = shape_create_deadstart(s)
        s['segs'] = self.get_segs()
        s['rcp'], s['zcp'] = self.seg_intersections(s['segs'], s['rb'], s['zb'])

        return s


    def update_plots(self, event=None):
        """
        METHOD: update_plots
        DESCRIPTION:                
        """      
                          
        try:                  
            s = self.update_data()

            # plot boundary shape
            for i in range(len(self.axs)):
                ax = self.axs[i]
                
                # remove old lines            
                if not self.hold_shape.get():
                    for artist in ax.lines + ax.collections + ax.texts:
                        artist.remove()
                
                # plot limiter and new shape
                self.plot_limiter(ax)
                ax.plot(s['rb'], s['zb'], linewidth=1, color='red')   
                
                # plot manually-defined (r,z) points
                for k in range(self.n_manual_control_pts):
                    rkey = 'r' + str(k+1)
                    zkey = 'z' + str(k+1)
                    ax.scatter(s[rkey], s[zkey], s=15, c='blue', alpha=1, marker='o')

                # plot x-points
                for k in range(self.n_xpts):
                    rkey = 'rx' + str(k+1)
                    zkey = 'zx' + str(k+1)
                    ax.scatter(s[rkey], s[zkey], s=50, c='red', alpha=1, marker='x')
                
                # plot control segments and points
                segs = s['segs']
                if self.plot_segs.get():
                    ax.scatter(s['rcp'], s['zcp'], s=15, c='blue', alpha=1, marker='.')
                    ax.plot(segs[:,[0,2]].T, segs[:,[1,3]].T, c='blue', alpha=0.3, linewidth=0.5)        
                
                # labels
                if self.label_control_pts.get():
                    for i in range(len(s['rcp'])):
                        txt = str(i+1)
                        ax.annotate(txt, (s['rcp'][i], s['zcp'][i]))
                
                if self.label_manual_control_pts.get():
                    for i in range(8):
                        txt = str(i+1)
                        rkey = 'r' + txt
                        zkey = 'z' + txt
                        ax.annotate(txt, (s[rkey], s[zkey]))

                if self.label_xpts.get():
                    for i in range(4):
                        txt = str(i+1)
                        rkey = 'rx' + txt
                        zkey = 'zx' + txt
                        ax.annotate(txt, (s[rkey], s[zkey]))

            self.canvas.draw()  

        except:
            # plot boundary shape
            for i in range(len(self.axs)):
                ax = self.axs[i]
                
                # remove old lines            
                if not self.hold_shape.get():
                    for artist in ax.lines + ax.collections + ax.texts:
                        artist.remove()
                
                # plot limiter
                self.plot_limiter(ax)
                
                        
        

    def load_shape_and_plot(self, filename=None, event=None):
        self.load_shape()
        self.update_plots()
   
    def load_shape_default(self, event=None):        
        self.load_shape(self.default_shape_fn)

    def load_shape_default_and_plot(self, event=None):
        self.load_shape_default()
        self.update_plots()

    def load_shape(self, filename=None, event=None):

        if filename is None:
            filetypes=[("JSON files","*.json"), ("Text Documents","*.txt"), ("All Files","*.*")]
            fid = tk.filedialog.askopenfile(filetypes=filetypes)
        else:
            fid = open(filename, 'r')

        s = json.load(fid)
        self.shape = s['shape_params']          
        
        # backwards compatibility
        if 'seg_params' in s.keys():
            for key in s['seg_params'].keys():
                self.shape[key] = s['seg_params'][key]

        self.shape = self.add_missing_shape_fields(self.shape)                
        self.shape = self.add_aux_geom_params(self.shape)

        for key in self.tk_shape.keys():
            self.tk_shape[key].set(str(self.shape[key]))                
                        
    def json_encode(self, d):
        for key in d:
            if isinstance(d[key], np.ndarray):
                d[key] = d[key].tolist()
            elif isinstance(d[key], dict):
                d[key] = self.json_encode(d[key])                                
        return d
    
    def save_file(self, d, filename=None):
        if filename is None:
            fid = tk.filedialog.asksaveasfile(initialfile='.json', defaultextension='.json',
                                        filetypes=[("All Files","*.*"),("Text Documents","*.txt")])                
        else:
            fid = open(filename, 'w')
        
        fid.write(json.dumps(d, indent=4))
        fid.close()

    def save_shape_default(self, event=None):                    
        self.save_shape(self.default_shape_fn)                

    def save_shape(self, filename=None, event=None):        
        s = self.update_data()          
        
        d = {'shape_params':s, 'geo':self.geo}
        d = self.json_encode(d)        
        self.save_file(d, filename)      
        print('Shape saved to file successfully.')  


def main():     

    parser = argparse.ArgumentParser(description='Specify a tokamak geometry to load')
    parser.add_argument('-tok', dest='tokamak', default='nstxu', 
        help="Specify a tokamak. Supported options are 'nstxu' [default]")
    args = parser.parse_args()

    app = ShapeApp(args.tokamak)
    app.root.mainloop()

if __name__ == '__main__':
    main()
