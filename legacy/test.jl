img = read(f)
open_video_out("/Users/cfranken/video3.mp4", Gray.(img[500:1201, 1550:2001]), framerate=framerate, encoder_options=encoder_options) do writer
    while !eof(f)
        img2 = read(f)
        write(writer, Gray.(img2[500:1201, 1550:2001])-Gray.(img[500:1201, 1550:2001]))
        img .= img2
    end
end