# tissue_classification

Code here first implement Bayesian classification with Gaussian smoothing in ITK, then extract the white matter(WM) from 
the segmentation. After the first step, the code run the marching cubes algorithm on your binary WM segmentation in VTK 
to obtain a surface representation of your segmentation (polydata). Finally, visualize the grayscale image and your polydata
with the vtkInteractorStyleTrackballCamera. 
