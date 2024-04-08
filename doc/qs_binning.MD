# Simulate quality score binning


You can simulate the base quality score binning by using the `--qs-bins <FILE>` option. The file should contain lines with the following format:
```
[RangeStart],[RangeEnd],[QualityScoreToAssign]
```
where the ranges are inclusive. All values should be comma-separated, >=0 and <=255. The ranges must cover the entire range from the first RangeStart to the last RangeEnd. The first range start must be 0, and the file is assumed to be sorted. 

Following is an example file for NovaSeq 6000 (RTA3) binning. Please note that these values may change depending on the sequencing platform and the provider.

```
0,2,2
3,14,12
15,30,23
31,40,37
```

Assigns quality score 2 to reads with quality scores 0-2, 12 to reads with quality scores 3-14, etc. Any quality score above 40 will be assigned 37. 

To assign a specific quality score to a single quality score value, simply use the same value for the range start and end. For example:

```
0,0,2
```

will assign quality score 2 to reads with quality score 0.