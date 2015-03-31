#!/bin/env perl

## Script to check two directory trees to ensure that they 
## are identical - this is useful to use after a complex merge
## if you want to ensure that the merged tree is identical
## to the development tree (e.g. after merging 'devel' into
## 'trunk' ensuring that 'trunk' is identical to 'devel')
##
## (C) Christopher Woods, 2006

sub compareDirs
{
      my ($dir0,$dir1) = @_;

      local(*DIR0);

      opendir DIR0,"$dir0" or die "Cannot read directory '$dir0': $!";

      my @files = readdir(DIR0);

      foreach my $file (@files)
      {
            if ($file =~ m/^\.|~$|\.o$/o)
            {
                  next;
            }

            my $fullfile0 = "$dir0/$file";
            my $fullfile1 = "$dir1/$file";
            
            if (not -e $fullfile1)
            {
                  print "$fullfile0 exists, but $fullfile1 doesn't!\n";
            }
            elsif ( -d $fullfile0 )
            {
                  compareDirs($fullfile0,$fullfile1);
            }
            else
            {
                  #diff the two files
                  $output = `diff $fullfile0 $fullfile1`;

                  if ($output)
                  {
                        print "$fullfile0 and $fullfile1 differ!\n";
                  }
            }
      }
      
      #finally, scan the contents of dir1 and ensure that they
      #all appear in dir0
      local(*DIR1);
      
      opendir DIR1,"$dir1" or die "Cannot open '$dir1': $!";
      
      @files = readdir(DIR1);
      
      foreach $file (@files)
      {
            if ($file =~ m/^\.|~$|\.o$/o)
            {
                  next;
            }

            $fullfile1 = "$dir1/$file";
            $fullfile0 = "$dir0/$file";
            
            if ( ! -e $fullfile0 )
            {
                  print "$fullfile1 exists but $fullfile0 doesn't!\n";
            }
      }
}


$dir0 = $ARGV[0];
$dir1 = $ARGV[1];

compareDirs($dir0,$dir1);
