# BezierGreen

## Introduction
This GitHub repository contains the source code for the implementation of the paper "Flexible 3D Cage-based Deformation via Green Coordinates on Bézier Patches". Our implementation builds upon the excellent survey work: [CageModeler](https://github.com/DanStroeter/CageModeler). We sincerely thank the authors for their contributions—our work truly stands on the shoulders of giants.

## Run
Our configuration environment is similar to [CageModeler](https://github.com/DanStroeter/CageModeler), using CMake to generate the project, which has cross-platform generation capabilities. Here, we specifically introduce the experience of configuring with CMake-gui, vcpkg, and Visual Studio in Windows.

First, clone the repository:
```bash
git clone https://github.com/Submanifold/BezierGreen.git
```

In the root directory of the project, run:
```bash
cmake-gui
```

Where is the source code: Choose the root directory of this project (where `CMakeLists.txt` and `vcpkg.json` are located).

Where to build the binaries: Create a `build` folder in the root directory and select this folder.

Then, click `Configure`, select `Visual Studio`. Click `Specify toolchain file for cross-compiling`, and specify the location of your vcpkg (you can download vcpkg to the project root directory). The required packages are in `vcpkg.json`, and using CMake-gui to select the toolchain should automatically download the required packages.  If you encounter errors related to Vulkan during the configuration process, you can download [Vulkan](https://vulkan.lunarg.com/sdk/home) from the following link and configure the system environment variables accordingly.



After that, click `Generate` in cmake-gui, and the Visual Studio project `cageDeformation3D.sln` will be generated in the `build` directory. Select the `RelWithDebInfo` mode, build the solution, and the `cageDeformation3D.exe` will be generated. If any DLL files are missing, simply place the corresponding files in the same directory as `cageDeformation3D.exe`.

The abbreviation of our method is BGC.

Run Command:
```bash
./cageDeformation3D.exe --model=/path/to/model.obj --cage=/path/to/source_bezier_cage.txt --cage-deformed=/path/to/target_bezier_cage.txt --BGC
```
For instance, the following two command lines:
```bash
./cageDeformation3D.exe --model=../../../data/ManHead/ManHead.obj --cage=../../../data/ManHead/ManHead_bezier_cage.txt --cage-deformed=../../../data/ManHead/ManHead_bezier_cage_deformed.txt --BGC
./cageDeformation3D.exe --model=../../../data/Botijo/Botijo.obj --cage=../../../data/Botijo/Botijo_bezier_cage.txt --cage-deformed=../../../data/Botijo/Botijo_bezier_cage_deformed.txt --BGC
```
are used to generate the teaser figure of our work.

The `data` folder contains some data, and the cage is currently stored in txt files, which are either 3rd-order Bézier triangles or tensor product Bézier patches. Therefore, each line contains 30 (10 control points x, y, z for Bézier triangle) or 48 (16 control points x, y, z for tensor product Bézier patch) values, arranged in the following order:

![Control Points](./control_points.png)

Additionally, if you want to use the cross-product Neumann term method mentioned in Appendix D, you can change `cross_product_BGC` to `true` in the main function (around line 127). We are optimizing the code, and trying to support inputs of different dimensions for better scalability, although using dynamic length for some variables is slightly slower than using fixed-length arrays for 3rd-order cases.

## Visualization
You can directly use [Meshlab](https://www.meshlab.net/) to visualize the generated obj file. To visualize the Bézier cage and the mesh together, we provide a [Blender](https://www.blender.org/) script named `vis.blend` in the `./vis` folder, which includes some code in the script. First, make sure that the `obj_filepath` and `bezier_filepath` on about line 270 are correct. Then click `Run script`, switch to the `Viewpoint Shading` mode, and you can view the result. You can adjust the style such as colors according to your preferences.

## Citation
```
@inproceedings{10.1145/3721238.3730630,
author = {Dong Xiao and Renjie Chen},
title = {Flexible 3D Cage-based Deformation via Green Coordinates on B\'{e}zier Patches},
year = {2025},
isbn = {9798400715402},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3721238.3730630},
doi = {10.1145/3721238.3730630},
series = {SIGGRAPH Conference Papers '25}
```
