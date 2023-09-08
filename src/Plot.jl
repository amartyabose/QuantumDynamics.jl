module Plot

using Reexport

@reexport using LaTeXStrings
import PyPlot
const plt = PyPlot

export plt, new_figure, contourf

function new_figure(plot_type="full")
    if plot_type == "full" || plot_type == "2D"
        fig = plt.figure(; figsize=(3.175, 2.25))
    elseif plot_type == "half"
        fig = plt.figure(; figsize=(1.56, 1.4))
    end
    ax = fig.add_subplot(111)

    for axis in ["top", "bottom", "left", "right"]
        ax.spines[axis].set_linewidth(0.1)
        ax.spines[axis].set_color("gray")
    end
    if plot_type == "2D"
        ax.tick_params(width=0.1, direction="out", color="gray")
    else
        ax.tick_params(width=0.1, direction="in", color="gray")
    end

    fig, ax
end

function contourf(x, y, z; levels=100, cmap="viridis")
    cnt = plt.contourf(x, y, z, levels)
    for c in cnt.collections
        c.set_edgecolor("face")
    end
    cbar = plt.colorbar(cnt, drawedges=false)
    cbar.outline.set_linewidth(0.1)
    cbar.outline.set_edgecolor("gray")
    cbar.ax.tick_params(width=0.1, direction="out", color="gray")
    cbar.solids.set_edgecolor("face")
    cnt, cbar
end

end