In recent years a multitude of parallel programming models have been introduced to ease parallel programming. Each programming model brings its own concepts and semantics, which makes it hard to see their impact on performance. Starbench is a benchmark suite that allows comparing different parallel programming models for embedded and consumer applications. Starbench consist of C/C++ benchmarks and currently covers video coding, image compression, image processing, hashing, artificial intelligence, computer vision, and compression. For each of the benchmark an optimized Pthreads version has been developed to serve as baseline. The suite has been succesfully used to evaluate the versatility and efficiency of OmpSs, a task-based programming model developed in the Encore project. While the main target is to provide a means to evaluate programming models and runtime improvements for E&C applications, Starbench and/or its individual benchmarks are also very suitable for other field of research such as computer architecture which require state-of-the-art parallel applications.

If you have any questions regarding Starbench, please write a mail to: starbench@aes.tu-berlin.de.

People involved:
* Michael Andersch
* Chi Ching Chi
* Prof. Dr. Ben Juurlink

All of the benchmarks included in UniBench were modified from its original source code, to allow us to compare perfomances with the serial version of the code.

Currently available benchmarks in UniBench are: kmeans, md5, ray-rot, rgbyuv, rotate, rotcc and streamcluster.

It is important to mention that the benchmarks were modified to use OpenMP task parallelism instead of OmpSs. And the programs were modified to run the serial version right after the parallel versio and display both times.

