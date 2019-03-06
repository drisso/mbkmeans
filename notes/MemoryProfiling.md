# Built-in R (base/utils package)

- Object size: `object.size` tells you the size of a single object. 
- `Rprof` captures usage frequently (Hadley: "profiling can, at best, capture memory usage every 1 ms"). Running `gctorture` forces gc after every memory allocation, which can help to slow down a program. 
- `tracemem` This function marks an object so that a message is printed whenever the internal code copies the object. `untracemem` undoes it.

# pryr

- `mem_used()` tells total size of all objects in memory. From Advanced R: 

>> This number won’t agree with the amount of memory reported by your operating system for a number of reasons:
>>
>> * It only includes objects created by R, not the R interpreter itself.
>> * Both R and the operating system are lazy: they won’t reclaim memory until it’s actually needed. R might be holding on to memory because the OS hasn’t yet asked for it back.
>> * R counts the memory occupied by objects but there may be gaps due to deleted objects. This problem is known as memory fragmentation.


- `mem_changed` tells you how memory changes during code execution, so that the input to the function would be an excecution.


# lineprof (Hadley)

This uses `Rprof`, but nicer feed back and pairs it up with your code. It requires that you use source to load the code, however. (http://adv-r.had.co.nz/memory.html#memory-profiling). 

>> Using lineprof is straightforward. source() the code, apply lineprof() to an expression, then use shine() to view the results. Note that you must use source() to load the code. This is because lineprof uses srcrefs to match up the code and run times. The needed srcrefs are only created when you load code from disk.

>> Next to the source code, four columns provide details about the performance of the code:
>>  * t, the time (in seconds) spent on that line of code (explained in measuring performance).
>>  * a, the memory (in megabytes) allocated by that line of code.
>>  * r, the memory (in megabytes) released by that line of code. While memory allocation is deterministic, memory release is stochastic: it depends on when the GC was run. This means that memory release only tells you that the memory released was no longer needed before this line.
>>  * d, the number of vector duplications that occurred. A vector duplication occurs when R copies a vector as a result of its copy on modify semantics.
>> You can hover over any of the bars to get the exact numbers. In this example, looking at the allocations tells us most of the story:

# profmem (Henrick)

# proftools (Tierney)

