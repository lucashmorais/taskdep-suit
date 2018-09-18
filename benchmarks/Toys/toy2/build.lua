project = {}

args = ""

project.build = function()
    os.execute(compiler .. " -O3 src/toy2.cpp -lm -fopenmp " .. flags .. " -o build/toy2-omp")
    os.execute(compiler .. " -O3 src/toy2.cpp -lm " .. flags .. " -o build/toy2-seq")
end

project.run = function()
    print("Sequential")
    os.execute("./build/toy2-seq " .. args .. " > output/seq");
    print("OpenMP")
    os.execute("./build/toy2-omp " .. args .. " > output/omp");
end

project.clean = function()
    os.execute("rm build/* output/*")
end

project[action]()