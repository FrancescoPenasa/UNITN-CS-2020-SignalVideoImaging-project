import os
import numpy as np

# take a nii file
data_path = "D:\dev\CT-Covid-19-August2020"
example_filename = os.path.join(data_path, 'volume-covid19-A-0000.nii.gz')

# load image
import nibabel as nib
img = nib.load(example_filename)

print ("img:")
print(img.shape)
