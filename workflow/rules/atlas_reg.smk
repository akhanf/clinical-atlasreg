
#if custom path is defined in the config for subject, use that
def get_t1w_filename(wildcards): 
    if wildcards.subject in config['subject_t1w_custom']:
        return config['subject_t1w_custom'][wildcards.subject]
    else:
        return config['subject_t1w']


rule import_subj_t1:
    input: get_t1w_filename
    output: bids(root='results',subject='{subject}',suffix='T1w.nii.gz')
    group: 'preproc'
    shell: 'cp {input} {output}'

rule rigid_aladin:
    input: 
        flo = bids(root='results',subject='{subject}',suffix='T1w.nii.gz'),
        ref = config['template_t1w'],
    output: 
        warped_subj = bids(root='results',subject='{subject}',suffix='T1w.nii.gz',space='{template}',desc='rigid'),
        xfm_ras = bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='rigid',type_='ras'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'reg_aladin -rigOnly -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.xfm_ras}'

rule affine_aladin:
    input: 
        flo = bids(root='results',subject='{subject}',suffix='T1w.nii.gz'),
        ref = config['template_t1w'],
    output: 
        warped_subj = bids(root='results',subject='{subject}',suffix='T1w.nii.gz',space='{template}',desc='affine'),
        xfm_ras = bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='affine',type_='ras'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.xfm_ras}'


rule qc_rigid:
    input:
        ref = config['template_t1w'],
        flo = bids(root='results',subject='{subject}',suffix='T1w.nii.gz',space='{template}',desc='rigid'),
    output:
        png = bids(root='qc',subject='{subject}',suffix='regqc.png',from_='subject', to='{template}',desc='rigid'),
        html = bids(root='qc',subject='{subject}',suffix='regqc.html',from_='subject', to='{template}', desc='rigid'),
    group: 'preproc'
    notebook: '../notebooks/vis_regqc.py.ipynb'

rule convert_xfm_ras2itk:
    input:
        bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='ras'),
    output:
        bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='itk'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'c3d_affine_tool {input}  -oitk {output}'

rule warp_mask_from_template:
    input: 
        mask = config['template_mask'],
        ref = bids(root='results',subject='{subject}',suffix='T1w.nii.gz'),
        xfm = bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='itk'),
    output:
        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',from_='{template}',reg='{desc}',desc='brain'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell: 'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} '
            ' -t [{input.xfm},1] ' #use inverse xfm (going from template to subject)

rule n4biasfield:
    input: 
        t1 = bids(root='results',subject='{subject}',suffix='T1w.nii.gz'),
        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',from_='MNI152NLin6Asym',reg='affine',desc='brain'),
    output:
        t1 = bids(root='results',subject='{subject}',desc='n4',suffix='T1w.nii.gz'),
    threads: 8
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}'


rule ants_syn_affine_init:
    input: 
        flo = bids(root='results',subject='{subject}',desc='n4',suffix='T1w.nii.gz'),
        ref = config['template_t1w'],
        flo_mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',from_='{template}',reg='affine',desc='brain'),
        ref_mask = config['template_mask'],
        init_xfm = bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='affine',type_='itk'),
    params:
        out_prefix = bids(root='results',suffix='_',from_='subject',to='{template}',subject='{subject}'),
        base_opts = '-d {dim} --float 1 --verbose 1 --random-seed {random_seed}'.format(dim=config['ants']['dim'],random_seed=config['ants']['random_seed']),
        intensity_opts = config['ants']['intensity_opts'],
        init_transform = lambda wildcards, input: '-r {xfm}'.format(xfm=input.init_xfm),
        masks = lambda wildcards, input: '--masks [{ref_mask},{flo_mask}]'.format(ref_mask=input.ref_mask,flo_mask=input.flo_mask),
        linear_multires = '-c [{reg_iterations},1e-6,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                reg_iterations = config['ants']['linear']['reg_iterations'],
                                shrink_factors = config['ants']['linear']['shrink_factors'],
                                smoothing_factors = config['ants']['linear']['smoothing_factors']),
        linear_metric = lambda wildcards, input: '-m MI[{template},{target},1,32,Regular,0.25]'.format( template=input.ref,target=input.flo),
        deform_model = '-t {deform_model}'.format(deform_model = config['ants']['deform']['transform_model']),
        deform_multires = '-c [{reg_iterations},1e-9,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                reg_iterations = config['ants']['deform']['reg_iterations'],
                                shrink_factors = config['ants']['deform']['shrink_factors'],
                                smoothing_factors = config['ants']['deform']['smoothing_factors']),
        deform_metric = lambda wildcards, input: '-m {metric}[{template},{target},1,4]'.format(
                                metric=config['ants']['deform']['sim_metric'],
                                template=input.ref, target=input.flo)

    output:
        out_affine = bids(root='results',suffix='_0GenericAffine.mat',from_='subject',to='{template}',subject='{subject}'),
        warped_flo = bids(root='results',suffix='T1w.nii.gz',space='{template}',desc='SyN',subject='{subject}'),
#    log: bids(root='logs',suffix='ants_SyN.log',subject='{subject}')
    threads: 16
    resources:
        mem_mb = 16000, # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time = 60 # 1 hrs
    container: config['singularity']['ants']
    group: 'preproc'
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsRegistration {params.base_opts} {params.intensity_opts} '
        '{params.init_transform} ' #initial xfm
        '-t Rigid[0.1] {params.linear_metric} {params.linear_multires} ' # rigid registration
        '-t Affine[0.1] {params.linear_metric} {params.linear_multires} ' # affine registration
        '{params.deform_model} {params.deform_metric} {params.deform_multires} '  # deformable registration
        '-o [{params.out_prefix},{output.warped_flo}]'

rule tissueseg_fsl:
    input:
        t1w = bids(root='results',subject='{subject}',desc='n4',suffix='T1w.nii.gz'),
    params:
        out_basename = bids(root='results',subject='{subject}',desc='n4',suffix='T1w_fast'),
    output:
        done = touch(bids(root='results',subject='{subject}',desc='n4',suffix='T1w_fast.done'))
    container: config['singularity']['neuroglia']
    group: 'tissueseg'
    shell:
        'fast -o {params.out_basename} --nobias --class=3 --verbose --segments {input}' 



#brainmasking with bet -- doesn't work well.. not using it..
rule bet:
    input:
        flo = bids(root='results',subject='{subject}',suffix='T1w.nii.gz')
    params:
        frac = 0.5,
        out_prefix = bids(root='results',subject='{subject}',from_='bet',suffix='desc-brain'),
    output: 
        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',desc='brain'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell: 'bet {input} {params.out_prefix} -f {params.frac} -m -n -v'


#rule qc_brainmask:
#    input: 
#        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',desc='brain'),
#        img = config['subject_t1w'],
#    output:
#        png = report(bids(root='qc',suffix='mask.png',desc='topup',**subj_wildcards),\
#            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
#
#        html = report(bids(root='qc',suffix='mask.html',desc='topup',**subj_wildcards),\
#            caption='../report/brainmask_dwi.rst', category='brainmask_dwi',\
#            subcategory=bids(**subj_wildcards,include_subject_dir=False,include_session_dir=False))
#    notebook: '../notebooks/vis_overlay_mask.py.ipynb'
        


#TODO ants, using aladin as initial rigid transform 



