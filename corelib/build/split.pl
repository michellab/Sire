#!/bin/env perl

while ($line = <STDIN>)
{
    @words = split(" ",$line);
    print "$words[$ARGV[0]]\n";
}
