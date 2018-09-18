project = {}

args = ""

q = "3"

project.build = function()
    os.execute(compiler .. " -O3 src/toy" .. q ..".cpp -lm -fopenmp " .. flags .. " -o build/toy" .. q .."-omp")
    os.execute(compiler .. " -O3 src/toy" .. q ..".cpp -lm " .. flags .. " -o build/toy" .. q .."-seq")
end

project.run = function()
    print("Sequential")
    os.execute("./build/toy" .. q .."-seq " .. args .. " > output/seq");
    print("OpenMP")
    os.execute("./build/toy" .. q .."-omp " .. args .. " > output/omp");
end

project.clean = function()
    os.execute("rm build/* output/*")
end

project[action]()