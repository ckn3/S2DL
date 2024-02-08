# S2DL: A Superpixel-based and Spatially-regularized Diffusion Learning Method for Unsupervised Hyperspectral Image Clustering

This code is an implementation of **Superpixel-based and Spatiallyregularized Diffusion Learning** proposed in "Superpixel-based and Spatially-regularized Diffusion Learning Method for Unsupervised Hyperspectral Image Clustering", see [Link](https://arxiv.org/abs/2312.15447). S2DL can be used as a clustering method for remote sensing datasets.

S2DL uses several Matlab Toolboxes, such as [Entropy Rate Superpixel](https://github.com/mingyuliutw/EntropyRateSuperpixel), and [Diffusion Learning](https://github.com/sampolk/DiffusionLearning).

Notes:
- Contact: kangnicui2@gmail.com

*If you find it useful or use it in any publications, please cite the following papers:*
## References
**Cui, K., Li, R., Polk, S.L., Lin, Y., Zhang, H., Murphy, J.M., & Plemmons, R. J., Chan, R. H.**. "Superpixel-based and Spatially-regularized Diffusion Learning Method for Unsupervised Hyperspectral Image Clustering". in *ArXiv*, 2023. [Link](https://arxiv.org/abs/2312.15447).

## Custom Data Integration Guide

To adapt the S2DL codebase for your own multispectral dataset, follow these steps:

### Step 1: Preparing Your Data

1. Access the data loading script at: [loadHSI.m](https://github.com/ckn3/S2DL/blob/main/backEnd/Misc/Preprocessing/loadHSI.m).
2. Add your dataset by inserting the following code within the `if-elseif` chain:

    ```matlab
    elseif strcmp(HSIName, 'Your_Data_Name')
        HSI = % Load your Multispectral data here, with shape m*n*p, where:
              % m and n are the spatial dimensions, and p denotes the number of spectral bands.
        GT = % Load your ground truth for the Multispectral data here, with shape m*n.
    ```

    Replace 'Your_Data_Name' with your dataset's identifier and provide the appropriate loading commands for `HSI` and `GT`.

### Step 2: Modifying the Main Script

1. Modify `S2DL_main.m` to include your dataset in the selection prompt. Update lines 9-17 as follows:

    ```matlab
    prompt = 'Which dataset? \n 1) Indian Pines (Corrected) \n 2) Salinas (Corrected) \n 3) Salinas A (Corrected) \n 4) Your Data \n';
    DataSelected = input(prompt);
    if DataSelected > 4 || DataSelected < 1
        disp('Incorrect prompt input. Please enter a valid number [1-4].')
    end

    datasets = {'IndianPines', 'Salinas', 'SalinasA', 'Your_Data_Name'};
    ```

    Ensure 'Your_Data_Name' matches the identifier used previously.

### Step 3: Running Your Dataset

Now, you're set to run S2DL with your dataset. Simply choose your data when prompted by `S2DL_main.m`.

By following these steps, you'll be able to integrate and utilize your own multispectral dataset within the S2DL framework, expanding its application to new and unique datasets.
