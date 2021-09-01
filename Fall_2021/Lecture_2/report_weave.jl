using Weave

filename = "/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_2/tennis_ball_hitting_door.jl"

# Julia markdown to HTML
weave(filename; doctype = "md2html", out_path = :pwd)
