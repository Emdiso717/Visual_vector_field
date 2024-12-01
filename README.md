# Visual vector field

## Compile

Compile with the following instructions：

    mkdir build
    cd build
    cmake ..
    make

## Run

Open the file `project.sln` to use the `Release` mode to run the program in `vs2022` or later.

## Parameter Adjustment

There is an `input.json` file in the project root file.Open the file to modify the program running parameters.

```json
  "path_of_mesh": "coherent1012/2.obj",
  "path_of_grad": "coherent1012/2.txt",

  "show_grad_line": 0,
  "grad_line_color": [ 0.0, 0.0, 0.0 ],
  "grad_line_length": 1,

  "point_color": [ 1.0, 0.0, 0.0 ],
  "num_points": 10,
  "point_size": 3,
  "length_of_point": 1,

  "moving_point_step": 0.5,
  "time_per_step": 0.01,
  "camera_zoom": 2,
  "time_restart": 10
```

- The `path of mesh` is the name of the `.obj` file that needs to be rendered.
- The `path of grad` is the gradient file to be read.

> [!TIP]
>
> **Place the files you want to read in under the `./data` folder in the root directory.This way you only need to fill in the name of the file.**

- The `show grad line` is used to control whether a vector direction line segment is displayed.
  - `1` is need to be displayed and `0` is not need.
- The `grad line color` indicates the color of the line segment.
- The `grad_line_length` indicates the length of the line segment.
- The `point_color` indicates the color of the moving points.
- The `num_points` indicates the number of the moving points.
- The `point_size` indicates the size of the moving points.
- The `length_of_point` indicates the length of the moving points.
- The `moving_point_step` indicates the step size of the moving points.
- The `time_per_step` control the time interval of each step.
- The `camera zoom` controls the zoom of the lens.
- The `time_restart` is the time to restart.

## File

The performance analysis and effects demonstration videos are stored in the `./doc&video` folder。
