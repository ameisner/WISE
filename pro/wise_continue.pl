#!/usr/bin/env perl
# 2012-Feb-29 - Written by Aaron Meisner

@proc = ("");

while (@proc) {
  @proc = `pgrep idl -u ameisner`;
  sleep(60);
}
