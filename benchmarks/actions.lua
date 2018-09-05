


build = function()
    os.execute("cd " ..  bench.group .. "/" .. bench.bench .. " && " .. lua .. " -e \'action = \"build\"; compiler = \"" .. compiler .. "\"; flags = \"" .. flags .. "\";\'" .. " build.lua")
end

run = function()
    os.execute("cd " ..  bench.group .. "/" .. bench.bench .. " && " .. lua .. " -e \'action = \"run\";\'" .. " build.lua")
end

clean = function()
    os.execute("cd " ..  bench.group .. "/" .. bench.bench .. " && " .. lua .. " -e \'action = \"clean\";\'" .. " build.lua")
end

flow = function(func)
    if arg[3] then
        bench.group = arg[2]
        bench.bench = arg[3]
        func()
    elseif arg[2] ~= "all" then
        for i, v in ipairs(benchmarks[arg[2]]) do
            bench.group = arg[2]
            bench.bench = v
            func()
        end
    else
        for k, t in pairs(benchmarks) do
            for i, v in ipairs(t) do
                bench.group = k
                bench.bench = v
                func()
            end
        end
    end
end

-- Actions implementations

buildBench = function()
    flow(build)
end

runBench = function()
    flow(run)
end

cleanBench = function()
    flow(clean)
end

listBenchs = function()
    print("List of available benchmarks:")
    for k, v in pairs(benchmarks) do
        for i, j in ipairs(v) do
            print("\t" .. k .. " " .. j)
        end
        print()
    end
end

helpMessage = function()
    print([[
bench is an alternative to UniBench
made by Vitor Silva, many thanks to UniBench developers

For now the implementation utilizes some globals, so please,
run the command on the root of the benchmarks folder.

./bench help ----------------- display this help message
./bench list ----------------- show all available benchmarks
./bench build all ------------ build all available benchmarks
./bench build group ---------- build all benchmarks from the group
./bench build group bench ---- build a specific benchmark
./bench run all -------------- run all available benchmarks
./bench run group ------------ run all benchmarks from the group
./bench run group bench ------ run a specific benchmark
./bench clean all ------------ clean all available benchmarks
./bench clean group ---------- clean all benchmarks from the group
./bench clean group bench ---- clean a specific benchmark

    ]])
end