
#include <libcpuid/libcpuid.h>

#include <stdio.h>

int main(int argc, char **argv)
{
    int i = 0;

    if (!cpuid_present())
    {
        printf("CPUID is not supported.\n");
        return -1;
    }
    
    struct cpu_raw_data_t raw_data;
       
    if (cpuid_get_raw_data(&raw_data) != 0)
    {
        printf("Could not get raw CPUID data.\n");
        return -1;
    }
      
    struct cpu_id_t cpuid;
        
    if (cpu_identify(&raw_data, &cpuid) != 0)
    {
        printf("Could not identify the CPU.\n");
        return -1;
    }
        
    printf("CPU information\n");
    printf("vendor = %s\n", cpuid.vendor_str);
    printf("brand = %s\n", cpuid.brand_str);
    printf("codename = %s\n", cpuid.cpu_codename);
    printf("num_cores = %d\n", cpuid.num_cores);
    printf("num_logical_cores = %d\n", cpuid.num_logical_cpus);
    printf("total_logical_cores = %d\n", cpuid.total_logical_cpus);
    printf("cpu_clock_os = %d MHz\n", cpu_clock_by_os());
    printf("cpu_clock = %d MHz\n", cpu_clock_measure(200,1));

    printf("\nCache information\n");
    printf("l1_data_cache = %d kilobytes\n", cpuid.l1_data_cache);
    printf("l1_instruction_cache = %d kilobytes\n", cpuid.l1_instruction_cache);
    printf("l2_cache = %d kilobytes\n", cpuid.l2_cache);
    printf("l3_cache = %d kilobytes\n", cpuid.l3_cache);

    printf("\nCPU supported features\n");

    for (i=0; i<NUM_CPU_FEATURES; ++i)
    {
        if (cpuid.flags[i])
            printf("%s ", cpu_feature_str(i));
    }

    printf("\n\nAll complete.\n");

    return 0;
}
