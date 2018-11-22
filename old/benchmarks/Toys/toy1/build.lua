project = {}

args = ""

project.build = function()
    os.execute(compiler .. " -O3 src/toy1.cpp -fopenmp " .. flags .. " -o build/toy1-omp")
    os.execute(compiler .. " -O3 src/toy1.cpp " .. flags .. " -o build/toy1-seq")
end

project.run = function()
    print("Sequential")
    os.execute("./build/toy1-seq " .. args .. " > output/seq");
    print("OpenMP")
    os.execute("./build/toy1-omp " .. args .. " > output/omp");
end

project.clean = function()
    os.execute("rm build/* output/*")
end

project[action]()