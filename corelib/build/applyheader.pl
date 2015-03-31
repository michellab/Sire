#!/usr/bin/perl -W

$tmpdir = "/tmp";

$templatefile = "build/templates/header";

sub readTemplate
{
    my ($file) = @_;
    
    local(*FILE);
    
    open FILE,"<$file" or die "Cannot read header $file: $!\n";
    
    my $i = 0;
    
    while (my $line = <FILE>)
    {
         $template[$i] = $line;
         $i++;
    }
    
    close(FILE);
}

sub getLastModifiedYear
{
      my ($filename) = @_;

      my @stats = stat($filename);
      my $modtime = $stats[9];
      my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($modtime);

      $year += 1900;
      return $year;
}

sub getAuthor
{
      my ($user) = @_;
      
      my %authors = ( "chris" => "Christopher Woods",
                      "chzcjw" => "Christopher Woods" );
                      
      if (defined $authors{$user})
      {
          return $authors{$user};
      }
      else
      {
          print "Who are you? :-> ";
          my $author = <STDIN>;
          chomp $author;
          
          print "Thanks - credit for the files will now be given to \"$author\".\n";
          
          $authors{$user} = $author;
          return $author;
      }
}

sub applyHeader
{
    my ($filename) = @_;
    
    local(*FILE);
    local(*NEWFILE);
    
    my $user = $ENV{USER};
    
    my $tmpfile = "$tmpdir/$user-applyheader.pl";

    open NEWFILE,">$tmpfile" or die "Cannot create $tmpfile: $!\n";

    open FILE,"<$filename" or die "Cannot open $file: $!\n";
    
    my $year = 0;
    my $author = 0;
    
    my $doneheader = 0;
    
    my $old_start = -1;
    my $old_end = -1;
    
    my $i = -1;
    
    #do a first read-through to get the existing header
    while (my $line = <FILE>)
    {
        $i++;
    
        if ($line =~ m/\/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\\/o)
        {
            #this is the first line of the header - keep reading in
            #the header to get the old year and author
            while (my $line = <FILE>)
            {
                $i++;
            
                chomp $line;
            
                if ($line =~ m/Copyright \(C\)\s+(\d\d\d\d)\s+([\w,\s]+)/o)
                {
                    $year = $1;
                    $author = $2;
                }
                elsif ($line =~ m/\\\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\//o)
                {
                    last;
                }
            }
        }
        elsif ($line =~ m/\@file/o)
        {
            #this contains an old-style header
            $old_start = $i-1;
            
            #read to the end of the old-style header
            while (my $line = <FILE>)
            {
                $i++;
                
                if ($line =~ m/\*\//o)
                {
                    $old_end = $i;
                    last;
                }
                elsif ($line =~ m/\(C\)\s+(\d\d\d\d)/o)
                {
                    $year = $1;
                }
            }
        }
    }
    
    close(FILE);
    
    #do we have an author and a year?
    if (not $year)
    {
        $year = getLastModifiedYear($filename);
    }
    
    if (not $author)
    {
        #look up the author from the username
        $author = getAuthor($user);
    }
    
    #ok - write the proper header to the tmp file
    if (scalar(@template) == 0)
    {
        readTemplate($templatefile);
    }
    
    foreach my $line (@template)
    {
        if ($line =~ m/Copyright \(C\)/o)
        {
            $line =~ s/<year>/$year/;
            $line =~ s/<name of author>/$author/;
        }
        
        print NEWFILE $line;
    }
    
    #now write the file to the tmp file (without its old header)
    open FILE,"<$filename" or die "Cannot open $file: $!\n";

    $i = -1;

    my $start = 1;

    while (my $line = <FILE>)
    {
        chomp($line);
    
        $i++;
    
        if ($start)
        {
            my @words = split(" ",$line);
            if (not defined $words[0])
            {
                next;
            }
        }
    
        if ($line =~ m/\/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\\/o)
        {
            #this is the first line of the header
            # -read away until we've reached the last line
            while (my $line = <FILE>)
            {
                if ($line =~ m/\\\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\//o)
                {
                    last;
                }
            }
        }
        elsif ($i == $old_start)
        {
            #read away the old style header!
            while (<FILE>)
            {
                $i++;
                
                if ($i == $old_end){ last;}
            }
        }
        else
        {
            $start = 0;
            
            print NEWFILE "$line\n";
        }
    }
    
    close(FILE);
    close(NEWFILE);
    
    #mv the tmp file over the original
    system("mv -f $tmpfile $filename");
}

foreach $file (@ARGV)
{
      applyHeader($file);
}
