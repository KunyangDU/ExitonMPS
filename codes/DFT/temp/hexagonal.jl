using CairoMakie

function hexagon_vertices(radius, center)
    angles = range(0, 2π, length=7)[1:6]
    vertices = Point2f[(center[1] + radius * cos(θ), center[2] + radius * sin(θ)) for θ in angles]
    return vertices
end

function plot_hexagon(radius, center)
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1], title = "Hexagon Plot")
    vertices = hexagon_vertices(radius, center)
    
    # 绘制六方形的多边形
    poly = poly!(ax, vertices, color = :blue)
    
    # 设置绘图区域的范围
    limits!(ax, center[1] - radius * 1.5, center[1] + radius * 1.5, center[2] - radius * 1.5, center[2] + radius * 1.5)
    
    # 显示绘图
    display(fig)
end

# 调用函数绘制六方形
plot_hexagon(1.0, (0.0, 0.0))
