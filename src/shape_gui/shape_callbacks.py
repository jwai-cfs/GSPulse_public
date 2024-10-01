import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from intersections import intersection

def shape_create_deadstart(s):
    """
    Create a shape based on the shaping parameters in s
    Parameters:
    - s (dict): A dictionary containing shaping parameters like triangularity, squareness, etc.
    Returns:
    - r (array): An array of radial coordinates.
    - z (array): An array of vertical coordinates.
    Description:
    This function creates a shape by generating a circle and adding x-points to it. 
    It then creates a convex hull from the circle and x-points, and interpolates to a higher point density. 
    The resulting shape is analyzed and shaped based on the given parameters. 
    Finally, the shape is interpolated and sorted to ensure a closed loop.
    """          

    # make a circle
    th = np.linspace(0, 2*np.pi, 400)
    x = np.cos(th)
    y = np.sin(th)

    # add the x-points
    x = np.append(x, 0)
    y = np.append(y, np.min(y) - s['c_xplo'])

    x = np.append(x, 0)
    y = np.append(y, np.max(y) + s['c_xpup'])


    # create a convex hull from the circle + x-points
    xy = np.vstack((x,y)).T
    hull = ConvexHull(xy)
    x = xy[hull.vertices,0]
    y = xy[hull.vertices,1]
    
    # interpolate to higher point density
    x, y = interparc(x, y, 500, mergeit=False, forceloop=False)
    x, y = sort_ccw(x, y)
    
    # apply shaping parameters
    [r, z] = shape_edit(x, y, s)
    s = shape_analysis(r,z)

    # interpolate and sort
    r, z = sort_ccw(r, z)
    r, z = interparc(r, z, 500, mergeit=False, forceloop=True)
    
    return r, z



def sort_ccw(x,y):
    """
    sort points in counter-clockwise order
    """
    cx = x.mean()
    cy = y.mean()    
    angles = np.arctan2(y-cy, -x+cx)
    i = np.argsort(-angles)
    return x[i], y[i]


def interparc(x, y, n=100, forceloop=False, mergeit=False):
    """
    Interpolate points along a curve such that the arc length between points is equal.
    """
    if forceloop:
        if (x[0] != x[-1]) or (y[0] != y[-1]):
            x = np.append(x, x[0])
            y = np.append(y, y[0])
                
    lens = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    lens = np.insert(lens, 0, 0)
    arclens = np.cumsum(lens)

    # evenly distributed arc lengths
    s = np.linspace(0, arclens[-1], n+1)

    # linearly interpolate points according to arclength
    x2 = np.zeros(n)
    y2 = np.zeros(n)

    for i in range(n):
        k = np.where(arclens > s[i])[0][0] - 1
        dk = (s[i] - arclens[k]) / (arclens[k+1] - arclens[k]) # remainder
        x2[i] = x[k] + dk * (x[k+1] - x[k])
        y2[i] = y[k] + dk * (y[k+1] - y[k])
        
    # merge original set of points with new set
    if mergeit:
        v1 = np.vstack((x2,y2,s[:-1]))
        v2 = np.vstack((x,y,arclens))
        v = np.hstack((v1,v2))        
        
        _, iuniq = np.unique(v[2,:], return_index=True)     # index for unique points
        isort = np.argsort(v[2,iuniq])   # sort according to arc length
        idx = iuniq[isort]   # index of sorted, unique points

        x2 = v[0,idx]
        y2 = v[1,idx]        

    # forceloop for the new set of points
    if forceloop:
        if (x2[0] != x2[-1]) or (y2[0] != y2[-1]):
            x2 = np.append(x2, x2[0])
            y2 = np.append(y2, y2[0])

    return x2, y2

def shape_analysis(r,z):
    """
    
    """
    r = np.asarray(r)
    z = np.asarray(z)
    [r,z] = sort_ccw(r,z)

    # find inner, outer, upper, lower points
    s = {}
    ii = np.argmin(r)
    s['ri'] = r[ii]
    s['zi'] = z[ii]

    io = np.argmax(r)
    s['ro'] = r[io]
    s['zo'] = z[io]

    iu = np.argmax(z)
    s['ru'] = r[iu]
    s['zu'] = z[iu]

    il = np.argmin(z)
    s['rl'] = r[il]
    s['zl'] = z[il]

    s['R0'] = (s['ro']+s['ri'])/2
    s['Z0'] = (s['zu']+s['zl'])/2
    s['aminor'] = (s['ro']-s['ri'])/2
    s['bminor'] = (s['zu']-s['zl'])/2
    s['k'] = s['bminor'] / s['aminor']
    s['epsilon'] = s['aminor'] / s['R0']
    s['aspectratio'] = 1/ s['epsilon']
    s['triu'] = (s['R0'] - s['ru']) / s['aminor']
    s['tril'] = (s['R0'] - s['rl']) / s['aminor']
    s['tri'] = (s['triu'] + s['tril']) / 2

    # order matters for the squareness inputs 
    # (outer/inner point should precede upper/lower point) 
    s['squo'] = squareness(s['ro'], s['zo'], s['ru'], s['zu'], r, z)
    s['sqlo'] = squareness(s['ro'], s['zo'], s['rl'], s['zl'], r, z)
    s['squi'] = squareness(s['ri'], s['zi'], s['ru'], s['zu'], r, z)
    s['sqli'] = squareness(s['ri'], s['zi'], s['rl'], s['zl'], r, z)

    return s


def squareness(r1, z1, r2, z2, r, z):
    """
    squareness definition from: 
    https://iopscience.iop.org/article/10.1088/0741-3335/55/9/095009/meta
    """
    A = r1 - r2
    B = z2 - z1
    tol = (r2-r1)*1e-8

    rellipse = np.linspace(r2-tol, r1+tol, 100)
    zellipse = z1 + np.sign(B)*np.sqrt(B**2 - ((B/A)*(rellipse - r2))**2)

    rseg = np.array([r2,r1])
    zseg = np.array([z1,z2])

    [rc, zc] = intersection(rseg, zseg, rellipse, zellipse)
    [rd, zd] = intersection(rseg, zseg, r, z)

    LOD = np.linalg.norm(np.array([rd - r2, zd - z1]))
    LOC = np.linalg.norm(np.array([rc - r2, zc - z1]))
    LCE = np.linalg.norm(np.array([rc - r1, zc - z2]))

    sq = (LOD - LOC) / LCE

    return sq


 
def shape_edit(r,z,s):
    """
    reference for squareness:
    https://iopscience.iop.org/article/10.1088/0741-3335/55/9/095009/meta

    The parameters in s that will edit the shape are: ['R0','Z0','a','k','triu',
    'tril','squo','squi','sqlo','sqli','xplo','xpup'}.  All other parameters are
    ignored during shape_edit. 
    """
    r = np.asarray(r)
    z = np.asarray(z)

    [r,z] = sort_ccw(r,z)
    s0 = shape_analysis(r,z)
    
    # shape edits from (R0, Z0)
    r = r + s['R0'] - s0['R0']
    z = z + s['Z0'] - s0['Z0']
  
    # shape edits from a
    r = s['R0'] + (r-s['R0']) * s['aminor'] / s0['aminor']
    
  
    # shape edits from k
    b0 = s0['aminor'] * s0['k']
    bminor = s['aminor'] * s['k']
    z = s['Z0'] + (z-s['Z0']) * bminor / b0
  
    # shape edits from (triu, tril)
    s0 = shape_analysis(r,z)
    
    ru = s['R0'] - s['aminor'] * s['triu']
    dru = ru - s0['ru']         # how much ru needs to move to match triu
    iu = np.where(z > s0['Z0'])[0]    
    f = np.interp(r[iu], [s0['ri'], s0['ru'], s0['ro']], [0,1,0])
    r[iu] = r[iu] + f * dru    # movement is propto dru and distance from ri,ro
  
    rl = s['R0'] - s['aminor'] * s['tril']
    drl = rl - s0['rl']        
    il = np.where(z < s0['Z0'])[0]    
    f = np.interp(r[il], [s0['ri'], s0['rl'], s0['ro']], [0,1,0])
    r[il] = r[il] + f * drl    
    
    # shape edits from squareness
    s0 = shape_analysis(r,z)

    # order matters for the edit_squareness inputs
    # (outer/inner point should precede upper/lower point)
    [r1,z1] = edit_squareness(s0['ro'], s0['zo'], s0['ru'], s0['zu'], s0['squo'], s['squo'], r, z)
    [r2,z2] = edit_squareness(s0['ri'], s0['zi'], s0['ru'], s0['zu'], s0['squi'], s['squi'], r, z)
    [r3,z3] = edit_squareness(s0['ro'], s0['zo'], s0['rl'], s0['zl'], s0['sqlo'], s['sqlo'], r, z)
    [r4,z4] = edit_squareness(s0['ri'], s0['zi'], s0['rl'], s0['zl'], s0['sqli'], s['sqli'], r, z)

    r = np.concatenate((r1,r2,r3,r4))
    z = np.concatenate((z1,z2,z3,z4))
    
    [r,z] = sort_ccw(r,z)
    
    # make it a loop
    r = np.append(r, r[0])
    z = np.append(z, z[0])

    return r, z


def edit_squareness(r1,z1,r2,z2,sqinput,sqtarget,r,z):

    bminor = z2 - z1
    a = r1 - r2

    # (x,y) is the (r,z) normalized to the quadrant 1 unit circle
    x = (r-r2) / a
    y = (z-z1) / bminor
    i = np.where( (x>=0) & (y>=0))[0]  # use only quadrant 1    
    x = x[i]
    y = y[i]
    [x,y] = sort_ccw(x,y)
    th = np.arctan2(x,y)

    # curveA: the normalized input curve
    curveA = {'x':x, 'y':y, 'th':th}
    
    # curveB: the superellipse that matches input squareness, see ref
    n = -np.log(2) / np.log(1/np.sqrt(2) + sqinput*(1-1/np.sqrt(2)))
    x = np.linspace(1,0,100)
    y =  (1 - x**n) ** (1/n)
    th = np.arctan2(y,x)
    curveB = {'x':x, 'y':y, 'th':th}

    # curveC: the superellipse that matches target squareness, see ref
    n = -np.log(2) / np.log(1/np.sqrt(2) + sqtarget*(1-1/np.sqrt(2)))
    x = np.linspace(1,0,100)
    y = (1 - x**n) ** (1/n)
    th = np.arctan2(y,x)
    curveC = {'x':x, 'y':y, 'th':th}

    # interpolate all curves according to angle
    curveB['x'] = np.interp(curveA['th'], curveB['th'], curveB['x'])
    curveB['y'] = np.interp(curveA['th'], curveB['th'], curveB['y'])
    curveC['x'] = np.interp(curveA['th'], curveC['th'], curveC['x'])
    curveC['y'] = np.interp(curveA['th'], curveC['th'], curveC['y'])

    # shift input curveA by the amount that the superellipse shifted
    dx = curveC['x'] - curveB['x']
    dy = curveC['y'] - curveB['y']
    
    # curveD: the normalized output curve
    curveD = {}
    curveD['x'] = curveA['x'] + dx
    curveD['y'] = curveA['y'] + dy

    # denormalize
    r = curveD['x']*a + r2
    z = curveD['y']*bminor + z1
    
    return r, z

