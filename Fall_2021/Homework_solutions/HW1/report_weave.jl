using Weave

filename = "/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Homework_solutions/HW1/HW1_ball_bouncing.jl"

# Julia markdown to HTML
weave(filename; doctype = "md2html", out_path = :pwd)
