#!perl -T

use 5.006;
use strict;
use warnings;
use Test::More;

BEGIN {
    eval { require Math::AnyNum };
    plan skip_all => "Math::AnyNum is not installed"
      if $@;
}

plan tests => 36;

use Math::MatrixLUP;
use Math::AnyNum qw(:overload);

{
#<<<
    my $A = Math::MatrixLUP->new([
        [1, 2],
        [3, 4],
    ]);

    my $B = Math::MatrixLUP->new([
        [-3, -8, 3],
        [-2,  1, 4],
    ]);

    is_deeply(($A * $B)->as_array, [
        [-7,  -6,  11],
        [-17, -20,  25],
    ]);
#>>>
}

{
#<<<
    my $A = Math::MatrixLUP->new([
          [1, 2],
          [3, 4],
          [5, 6],
          [7, 8]
    ]);

    my $B = Math::MatrixLUP->new([
          [1, 2, 3],
          [4, 5, 6]
    ]);

    is_deeply($A->mul($B)->as_array, [
        [ 9,  12,  15],
        [19,  26,  33],
        [29,  40,  51],
        [39,  54,  69],
    ]);
#>>>
}

{
#<<<

    my $A = Math::MatrixLUP->new([
        [1,  1,  1,   1],
        [2,  4,  8,  16],
        [3,  9, 27,  81],
        [4, 16, 64, 256],
    ]);

    my $B = Math::MatrixLUP->new([
        [  4  , -3  ,  4/3,  -1/4 ],
        [-13/3, 19/4, -7/3,  11/24],
        [  3/2, -2  ,  7/6,  -1/4 ],
        [ -1/6,  1/4, -1/6,   1/24],
    ]);

    is_deeply(($A*$B)->as_array, Math::MatrixLUP->I(4)->as_array);
#>>>
}

{
#<<<
    my $A = Math::MatrixLUP->new([
        [2, 9, 4],
        [7, 5, 3],
        [6, 1, 8],
    ]);
#>>>

    is($A->det, -360);
}

{
#<<<
    my $A = Math::MatrixLUP->new([
        [  0,  1,  2,  3,  4 ],
        [  5,  6,  7,  8,  9 ],
        [ 10, 11, 12, 13, 14 ],
        [ 15, 16, 17, 18, 19 ],
        [ 20, 21, 22, 23, 24 ]
    ]);

    is($A->det, 0);
    is_deeply(($A**(-1))->as_array, $A->inv->as_array);
#>>>
}

{
#<<<
    my $A = Math::MatrixLUP->new([
        [1, 2, 0],
        [0, 3, 1],
        [1, 0, 0],
    ]);

    is_deeply($A->pow(8)->as_array, [
        [1291,  9580,  2930],
        [1465, 10871,  3325],
        [395,  2930,   896],
    ]);

    is_deeply($A->pow(9)->as_array, [
        [4221, 31322,  9580],
        [4790, 35543, 10871],
        [1291,  9580,  2930],
    ]);

    is_deeply(($A**10)->as_array, [
        [13801, 102408,  31322],
        [15661, 116209,  35543],
        [4221,  31322,   9580],
    ]);

    is_deeply(($A**-10)->as_array, $A->inv->pow(10)->as_array);
    is_deeply(($A**-10)->as_array, $A->pow(10)->inv->as_array);
#>>>
}

{
#<<<
    my $A = Math::MatrixLUP->new([
        [1, 1],
        [1, 0]
    ]);

    my $B = $A**100;
    is($B->[1][1] + $B->[1][0], $B->[0][0]);
    is($B->[0][0] - $B->[0][1], $B->[1][1]);
    is($B->[0][1], 354224848179261915075);
#>>>
}

{
    my $A = Math::MatrixLUP->new([[3, 1, 4], [1, 5, 9]]);
    my $B = Math::MatrixLUP->new([[2, 7, 1], [8, 2, 2]]);

    is_deeply(($A + $B)->as_array, [[5, 8,  5], [9,  7, 11]]);
    is_deeply(($A - $B)->as_array, [[1, -6, 3], [-7, 3, 7]]);

    is_deeply(($A + 42)->as_array, [[45,     43,     46],     [43,     47,     51]]);
    is_deeply(($A - 42)->as_array, [[-39,    -41,    -38],    [-41,    -37,    -33]]);
    is_deeply(($A * 42)->as_array, [[126,    42,     168],    [42,     210,    378]]);
    is_deeply(($A / 42)->as_array, [[1 / 14, 1 / 42, 2 / 21], [1 / 42, 5 / 42, 3 / 14]]);
    is_deeply(($A % 3)->as_array,  [[0,      1,      1],      [1,      2,      0]]);

    is_deeply(($A & 12)->as_array, [[0,  0,  4],  [0,  4,  8]]);
    is_deeply(($A | 3)->as_array,  [[3,  3,  7],  [3,  7,  11]]);
    is_deeply(($A ^ 42)->as_array, [[41, 43, 46], [43, 47, 35]]);
    is_deeply(($A << 3)->as_array, [[24, 8,  32], [8,  40, 72]]);
    is_deeply(($A >> 2)->as_array, [[0,  0,  1],  [0,  1,  2]]);

#<<<
    is_deeply($A->map(sub{my ($i, $j) = @_; $_ * $B->[$i][$j] })->as_array, [[6, 7, 4], [8, 10, 18]]);
    is_deeply($A->map(sub{my ($i, $j) = @_; $_ / $B->[$i][$j] })->as_array, [[3/2, 1/7, 4], [1/8, 5/2, 9/2]]);
    is_deeply($A->map(sub{my ($i, $j) = @_; $_ ** $B->[$i][$j] })->as_array, [[9, 1, 4], [1, 25, 81]]);
#>>>
}

{
#<<<
    my $A = Math::MatrixLUP->new([
        [2, -1,  5,  1],
        [3,  2,  2, -6],
        [1,  3,  3, -1],
        [5, -2, -3,  3],
    ]);

    my $B = Math::MatrixLUP->new([
        [1,  1,  1,   1],
        [2,  4,  8,  16],
        [3,  9, 27,  81],
        [4, 16, 64, 256],
    ]);

    is_deeply(($A/$B)->as_array, ($A*$B->inv)->as_array);
    is_deeply(($B/$A)->as_array, ($B*$A->inv)->as_array);

    is_deeply(('42' / $A)->as_array, ($A->inv*42)->as_array);
    is_deeply(('13' / $B)->as_array, ('13'*('1'/$B))->as_array);

    is_deeply(('42' - $A)->as_array, (-$A + 42)->as_array);
    is_deeply(('13' - $B)->as_array, (-$B + 13)->as_array);
    is_deeply(('-42' - $A)->as_array, ($A->neg - 42)->as_array);
#>>>
}
