clear all; clc; close all

shot = 204660;  
t = 0.1;

ROOT = getenv('RAMPROOT');
tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';
opts.cache_dir = [ROOT '/fetch/cache/'];
opts.plotit = 1;
eqs = fetch_eqs_nstxu(shot, t, tree, tokamak, server, opts);
eq = eqs.gdata;



