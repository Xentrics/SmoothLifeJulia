using Images, IJulia, Plots
gr() # FASTER!


grad = ColorGradient(:grays)
# grad = ColorGradient(RGBA{Float64}[RGBA{Float64}(0.05,0.05,0.05,1.0), RGBA{Float64}(0.95,0.95,0.95,1.0)], [0.0, 1.0])
heatmap(res["last"],
    c = grad, # seems to be the same as :grays, tbh.
    legend = nothing,
    xticks = nothing,
    yticks = nothing,
    aspect_ratio = 1.0)


function test_anim4()
    anim = @animate for i in 1:150
        plot(1:i, 1:i, color_palette = :grays)
    end

    gif(anim, "line_anim.gif", fps = 15)
    gif(anim, "line_anim.mp4", fps = 15)
    mp4(anim, "line_anim_.mp4", fps = 15)
    gif(anim, "line_anim.mov", fps = 15)
end


"""
Julia animate version
"""
function test_anim(curGrid)
    # define correct backend for plotting
    pyplot()
    # create animation object that can be saved as mp4 or displayed as gif
    anim = @animate for i=1:20
        #plt[:figure](); plt[:axis]("off"); plt[:imshow](curGrid, cmap=cm.Greys_r, vmin = 0.0, vmax = 1.0)
        Plots.plot(rand(5,5),linewidth=2,title="My Plot")
    end

    gif(anim, "test_anim.rand.gif", fps = 3)

    anim = @animate for i=1:20
        #plt[:figure](); plt[:axis]("off"); plt[:imshow](curGrid, cmap=cm.Greys_r, vmin = 0.0, vmax = 1.0)
        Plots.plot(curGrid, linewidth=2, title="My Plot")
    end

    gif(anim, "test_anim.curGrid.gif", fps = 3)
    return(anim);
end



function test_anim2(curGrid)
    function init()
    end

    function animate(i)
        im[:set_data](curGrid)
        plt[:draw]()
    end

    # create base figure
    fig = figure("MyFigure", figsize=(512, 512))

    myanim = anim.FuncAnimation(im, animate, init_func = init, frames = 100, interval=20)
    return(myanim)
end



# Does not seem to work for the real deal
function test_anim3(curGrid; nframes = 20, save_type = "mp4")
    function showanim(filename)
        base64_video = base64encode(open(filename))
        display("text/html", """<video controls src="data:video/x-m4v;base64,$base64_video">""")
    end

    # set up writer for movie / image files
    if save_type == "mp4"
        f_ending = ".mp4"
        #writer = anim.FFMpegWriter(fps=15, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
        writer = anim.FFMpegWriter(fps=15)
    elseif save_type == "gif"
        f_ending = ".gif"
        #writer = "imagemagick"
        writer = anim.ImageMagickWriter()
    else
        throw(error("Wrong save type!"))
    end



    fig = figure(figsize=(4,4))
    axis("off")

    function make_frame(i)
        #curGrid = min.(1, curGrid .* i)
        imshow(curGrid, cmap=cm.Greys_r, vmin = 0.0, vmax = 1.0)
    end

    withfig(fig) do
        myanim = anim.FuncAnimation(fig, make_frame, frames=nframes, interval=20)
        myanim[:save]("test3" * f_ending, dpi = 100, writer=writer)
    end
    close()
end
