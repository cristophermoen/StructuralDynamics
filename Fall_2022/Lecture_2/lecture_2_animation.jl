using Luxor, Colors

demo = Movie(400, 400, "test")

function backdrop(scene, framenumber)
    background("black")
end

function frame(scene, framenumber)
    sethue(Colors.HSV(framenumber, 1, 1))
    eased_n = scene.easingfunction(framenumber, 0, 1, scene.framerange.stop)
    circle(polar(100, -π/2 - (eased_n * 2π)), 80, :fill)
    text(string("frame $framenumber of $(scene.framerange.stop)"),
        Point(O.x, O.y-190),
        halign=:center)
    text(scene.opts,
        boxbottomcenter(BoundingBox()),
        halign=:center,
        valign=:bottom)
end

animate(demo, [
    Scene(demo, backdrop, 0:359),
    Scene(demo, frame, 0:359,
        easingfunction=easeinoutcubic,
        optarg="made with Julia")
    ],
    creategif=true)



# img = readpng("/Users/crismoen/Documents/dev/StructuralDynamics/Fall_2022/Lecture_2/137-1375722_carrinhos-hot-wheels-png-carros-da-hot-wheels.png")
# w = img.width
# h = img.height
# rulers()
# scale(0.3, 0.3)
# rotate(π/4)
# placeimage(img, Point(-w/2, -h/2), .5)


Drawing(300mm, 300mm, "post.svg")
	origin()
background("white")
	setline(2)
	sethue("grey")
box(O, 50mm, 10mm, action = :fill)
	sethue("white")
finish()
preview()	