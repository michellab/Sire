#!/usr/bin/perl -W

foreach $file (@ARGV)
{
    system("touch $file");
    system("perl ./build/applyheader.pl $file");
}

