from nilearn import plotting
import matplotlib.pyplot as plt


import matplotlib
matplotlib.use('Agg')

html_view = plotting.view_img(stat_map_img=snakemake.input.seg,bg_img=snakemake.input.img,
                              opacity=0.5,cmap='viridis',dim=-1,threshold=0.5,
                              symmetric_cmap=False,title='sub-{subject}'.format(**snakemake.wildcards))

html_view.save_as_html(snakemake.output.html)



display = plotting.plot_anat(snakemake.input.img,display_mode='ortho')
display.add_overlay(snakemake.input.seg,threshold=0.5)
display.savefig(snakemake.output.png)
display.close()

