#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 14:21:05 2024

@author: sei029
"""


def plot(project, site_ids, plotConvergence = True):
    
    # project = 'example'
    
    folder = f'../projects/{project}'
    
    import os
    
    # # extract site names
    # site_ids = []
    # for file in os.listdir(f'{folder}/transdMT/outfolder/csv'):
    #     if file.endswith(".csv") and not file.endswith("log.csv"):
    #         site_ids.append(file[:-4])
    
    
    def combine(site_id, folder):
        
        from PIL import Image
        #Read images
        image_fit = Image.open(f'{folder}/transdMT/outfolder/plots/fits/{site_id}_fit.png')
        image_model = Image.open(f'{folder}/transdMT/outfolder/plots/models/{site_id}_model.png')
        if plotConvergence:
            image_conver = Image.open(f'{folder}/transdMT/outfolder/plots/convergence/{site_id}_convergence.png')
        image_stats = Image.open(f'{folder}/transdMT/outfolder/plots/stats/{site_id}_stats.png')
        image_dim = Image.open(f'{folder}/edi/{site_id}.png')
        # image_dim.show()
        
        # #resize
        image_model_size = image_model.size
        image_fit_size = image_fit.size
        if plotConvergence:
            image_conver_size = image_conver.size
        image_stats_size = image_stats.size
        image_dim_size = image_dim.size
        
        if plotConvergence:
            image_conver = image_conver.resize((image_model_size[0], int(image_conver_size[1] * (image_model_size[0]/image_conver_size[0]))))
            image_conver_size = image_conver.size
        else:
            image_conver_size = 1,1,
        image_fit = image_fit.resize((int(image_fit_size[0]*(0.85*image_model_size[1]/image_fit_size[1])),int(0.85*image_model_size[1])))
        image_fit_size = image_fit.size
        image_stats = image_stats.resize((image_model_size[0], int(image_stats_size[1] * (image_model_size[0]/image_stats_size[0]))))
        image_stats_size = image_stats.size
        image_dim = image_dim.resize((int(image_dim.size[0]/2), int(image_dim.size[1]/2)))
        image_dim_size = image_dim.size
        
        # combine
        new_image = Image.new('RGB',(int(2.7*image_model_size[0]), image_model_size[1] + image_conver_size[1]), (250,250,250))
        new_image.paste(image_model,(0,0))
        if plotConvergence:
            new_image.paste(image_conver,(0,image_model_size[1]))
        new_image.paste(image_fit,(int(image_model_size[0]+0.06*image_model_size[0]),0))
        new_image.paste(image_stats,(image_conver_size[0],image_fit_size[1] + (image_model_size[1] + image_conver_size[1] - image_fit_size[1] - image_stats_size[1])))
        new_image.paste(image_dim, (int(image_model_size[0]+0.06*image_model_size[0]+image_fit_size[0]+100),100))
        new_image.save(f'{folder}/transdMT/outfolder/plots/combined/{site_id}.png',"PNG")
        # new_image.show()
        
    for site_id in site_ids:
        for filename in os.listdir(f'{folder}/transdMT/outfolder/csv'):
            if site_id in filename:
                print(f' Combining plots for MT site {site_id}...')
                combine(site_id, folder)
