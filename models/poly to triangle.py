input = "Sphere/abstract_ball.obj"
output = "Sphere/ball.obj"

with open(input, "r") as fileI:
    with open(output, "w") as fileO:
        for line in fileI:
            if line[0] != "f":
                fileO.write(line)
            else:
                if len(line.split()) != 5:
                    fileO.write(line)
                else:
                    caca = line.split()[1:]
                    fileO.write(f"f {caca[0]} {caca[1]} {caca[2]}\n")
                    fileO.write(f"f {caca[0]} {caca[2]} {caca[3]}\n")