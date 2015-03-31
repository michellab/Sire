#!/usr/bin/perl -W

foreach $file (@ARGV)
{
    system("touch $file");
    system("./build/applyheader.pl $file");
}

