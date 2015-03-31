
opendir DIR,".";

if (defined $ARGV[0])
{
    $spaces = $ARGV[0];
}
else
{
    $spaces = 6;
}

$indent = "";

for ($i=0; $i<$spaces; ++$i)
{
    $indent .= " ";
}

foreach $file (sort readdir(DIR))
{
    if ($file =~ m/\.cpp$/io)
    {
        print "$indent$file\n";
    }
}
