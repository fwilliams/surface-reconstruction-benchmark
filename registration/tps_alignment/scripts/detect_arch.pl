#!/usr/bin/perl -w

# Copyright 2007 Benedict Brown
#
# This script detects the platform on which the system is running, and
# sets appropriate compiler flags.  gcc should be getting an
# optimization option soon to do this, at which point we should move
# back to the trimesh2 system for detecting the architecture.

# this is to deal with Mac OS
$_ = `uname -a`;
if (/Darwin/) {
  
}

$_ = `cat /proc/cpuinfo`;

/vendor_id.*:.(.*)/;
$vendor = $1;

# /model.*:.(.*)/;
# $model = $1;

/model name.*:.(.*)/;
$name = $1;

# /stepping.*:.(.*)/;
# $stepping = $1;

/flags.*: (.*)/;
$flags = $1;

$arch = `uname -m`;

if ($arch =~ /x86_64/) {
    if ($vendor eq "AuthenticAMD") {
        print "opteron";
    } else { # intel
        print "nocona";
    }
} elsif ($arch =~ /i?86/) {
    if ($name =~ /Pentium.R. M/) {
        print "pentium-m";
    } elsif ($name =~ /Pentium.R. 4/ || $name =~ /Intel.R. XEON/) {
        if ($flags =~ /sse2/) {
            print "sse2";
        } elsif ($flags =~ /sse/) {
            print "pentium4";
        }
    }
}
