#Taskdep-suit

Taskdep-suit is a collection of benchmarks to test task dependecy based parallel programming models.

At the moment this suit aims to provide support for OMPSs.
Future support for OMPSs 2 and OMP4 is planned.

It may be used with ompss-testbench, by cloning it into the src folder: $/src/taskdep-suit

The conf file it the root of this project should make it easy to tune the compilation process, but is not used in all targets at the moment.

It requires the bench tool written in rust in order to provide the correct execution of the compilation process, but manual execution of makefiles may do the job.

State of ported benchmarks:

Kastors:
    - Jacobi:
        Working
    - SparseLu:
        Working, be carefull with values, it a demanding benchmark
    - Strassen:
        Working, has nested-tasks :(

Parsec:


TODO:

bench API needs to update the timing function to support RISC-V.