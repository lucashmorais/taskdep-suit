project = {}

args = "10000 50"

project.build = function()
    os.execute(compiler .. " -std=c99 -O3 src/pi.c src/test.c -lm -fopenmp " .. flags .. " -o build/pi-calc-omp")
    os.execute(compiler .. " -std=c99 -O3 src/pi.c src/test.c -lm " .. flags .. " -o build/pi-calc-seq")
end

project.run = function()
    print("Sequential")
    os.execute("./build/pi-calc-seq " .. args .. " > output/seq");
    print("OpenMP")
    os.execute("./build/pi-calc-omp " .. args .. " > output/omp");
end

project.clean = function()
    os.execute("rm build/* output/*")
end

project[action]()