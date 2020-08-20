# brainvisa_pipelines

The files here are used in multiple different pipelines.



Extracting sulci from BrainVISA is easiest with the UI, which takes ~6 minutes per brain and can run up to 20 brains in parallel. However, custom pipelines 
had to be built to handle higher resolution brains (i.e. greater resolution than 1mm), and have much longer processing times (i.e. 1 hour per hemisphere).
Parallelizing this code (or simply running it at different starting points in multiple terminals) will reduce run time.

To run sulci extraction, sulci labeling, and OFC sulcus isolation on 0.7mm brains, first enter a directory in your terminal with
a directory setup as subj>subj> (T1w.nii.gz + ribbon.nii.gz). In the top of this directory, type run_brainvisa.sh. This will begin iterating over 
all subdirectories to search for T1w.nii.gz and ribbon.nii.gz files to begin importing freesurfer segmentations to BrainVISA. Make sure you have BrainVISA installed
in your home directory, as well as freesurfer installed in /usr/local with the proper ~/.cshrc lines pointing to it (see freesurfer's documentation for install).
Additionally, make sure the brainvisa_files folder is in the proper directory, as these many files are accessed to apply to the pipeline at different steps.


After the pipeline, to manually extract the LOS and MOS, you must edit the names of sulci in anatomist (see video).
Then, run the script to extract the sulci to a version visible in itksnap.


The itksnap_img.nii files can be directed to isomap projections in the sulci_reduction code. This code helps present the moving averages along the projected
first axis.
