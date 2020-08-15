from nilearn import plotting
import matplotlib.pyplot as plt


import matplotlib
matplotlib.use('Agg')

html_view = plotting.view_img(stat_map_img=snakemake.input.seg,bg_img=snakemake.input.img,
                              opacity=0.5,cmap='viridis',dim=-1,threshold=0.1,
                              symmetric_cmap=False,title='sub-{subject}'.format(**snakemake.wildcards))

html_view.save_as_html(snakemake.output.html)



fig = plt.figure(figsize=(30, 3), facecolor='k')
display = plotting.plot_anat(snakemake.input.img,figure=fig,dim=0,display_mode='y')
display.add_overlay(snakemake.input.seg)
display.savefig(snakemake.output.png)
display.close()

