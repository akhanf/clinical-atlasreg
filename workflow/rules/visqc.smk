
rule qc_reg:
    input:
        ref = config['template_t1w'],
        flo = bids(root='results',subject='{subject}',suffix='T1w.nii.gz',space='{template}',desc='{desc}'),
    output:
        png = report(bids(root='qc',subject='{subject}',suffix='regqc.png',from_='subject', to='{template}',desc='{desc}'),
                caption='../reports/regqc.rst',
                category='Registration QC',
                subcategory='{desc} {template}'),
        html = report(bids(root='qc',subject='{subject}',suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
                caption='../reports/regqc.rst',
                category='Registration QC',
                subcategory='{desc} {template}'),
    group: 'preproc'
    script: '../scripts/vis_regqc.py'


rule qc_probseg:
    input:
        img = bids(root='results',subject='{subject}',desc='n4',from_='{template}', suffix='T1w.nii.gz'),
        seg = bids(root='results',subject='{subject}',suffix='probseg.nii.gz',label='{tissue}',from_='{template}'),
    output:
        png = report(bids(root='qc',subject='{subject}',suffix='probseg.png',label='{tissue}', from_='{template}'),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='Probabilistic Seg {tissue} from {template}'),
        html = report(bids(root='qc',subject='{subject}',suffix='probseg.html',label='{tissue}', from_='{template}'),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='Probabilistic Seg {tissue} from {template}'),
    group: 'preproc'
    script: '../scripts/vis_segqc.py'

rule qc_dseg:
    input:
        img = bids(root='results',subject='{subject}',desc='n4',from_='{template}', suffix='T1w.nii.gz'),
        seg = bids(root='results',subject='{subject}',suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}'),
    output:
        png = report(bids(root='qc',subject='{subject}',suffix='dseg.png',atlas='{atlas}', from_='{template}'),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='{atlas} Atlas from {template}'),
        html = report(bids(root='qc',subject='{subject}',suffix='dseg.html',atlas='{atlas}', from_='{template}'),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='{atlas} Atlas from {template}'),
    group: 'preproc'
    script: '../scripts/vis_segqc.py'



