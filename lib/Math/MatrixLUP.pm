package Math::MatrixLUP;

use 5.010;
use strict;
use warnings;

our $VERSION = '0.01';

use overload
  '""'  => \&stringify,
  'neg' => \&neg,

  '==' => \&eq,
  '!=' => \&ne,

  '+' => \&add,
  '*' => \&mul,

  '/' => sub { @_ = ($_[1], $_[0]) if $_[2]; goto &div },
  '-' => sub { @_ = ($_[1], $_[0]) if $_[2]; goto &sub },

  '**' => sub { @_ = $_[2] ? @_[1, 0] : @_[0, 1]; goto &pow },
  '%'  => sub { @_ = $_[2] ? @_[1, 0] : @_[0, 1]; goto &mod },
  ;

sub new {
    my ($class, $matrix) = @_;
    bless {A => $matrix}, $class;
}

sub _LUP_decomposition {
    my ($self) = @_;

    my @A = map { [@$_] } @{$self->{A}};
    my $N = $#A;
    my @P = (0 .. $N + 1);

    foreach my $i (0 .. $N) {

        my $maxA = 0;
        my $imax = $i;

        foreach my $k ($i .. $N) {
            my $absA = CORE::abs($A[$k][$i] // return [$N, \@A, \@P]);

            if ($absA > $maxA) {
                $maxA = $absA;
                $imax = $k;
            }
        }

        if ($imax != $i) {

            @P[$i, $imax] = @P[$imax, $i];
            @A[$i, $imax] = @A[$imax, $i];

            ++$P[$N + 1];
        }

        foreach my $j ($i + 1 .. $N) {

            if ($A[$i][$i] == 0) {
                return [$N, \@A, \@P];
            }

            $A[$j][$i] /= $A[$i][$i];

            foreach my $k ($i + 1 .. $N) {
                $A[$j][$k] -= $A[$j][$i] * $A[$i][$k];
            }
        }
    }

    [$N, \@A, \@P];
}

sub decompose {
    my ($self) = @_;
    $self->{_decomposed} //= $self->_LUP_decomposition;
}

sub clone {
    my ($self) = @_;
    __PACKAGE__->new([map { [@$_] } @{$self->{A}}]);
}

sub transpose {
    my ($self) = @_;

    my $A = $self->{A};

    my $rows = $#{$A};
    my $cols = $#{$A->[0]};

    __PACKAGE__->new(
        [
         map {
             my $i = $_;
             [map { $A->[$_][$i] } 0 .. $rows]
           } 0 .. $cols
        ]
    );
}

sub scalar_mul {
    my ($matrix, $scalar) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        push @B, [map { $_ * $scalar } @$row];
    }

    __PACKAGE__->new(\@B);
}

sub scalar_add {
    my ($matrix, $scalar) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        push @B, [map { $_ + $scalar } @$row];
    }

    __PACKAGE__->new(\@B);
}

sub scalar_sub {
    my ($matrix, $scalar) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        push @B, [map { $_ - $scalar } @$row];
    }

    __PACKAGE__->new(\@B);
}

sub scalar_div {
    my ($matrix, $scalar) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        push @B, [map { $_ / $scalar } @$row];
    }

    __PACKAGE__->new(\@B);
}

sub scalar_mod {
    my ($matrix, $scalar) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        push @B, [map { $_ % $scalar } @$row];
    }

    __PACKAGE__->new(\@B);
}

sub neg {
    my ($matrix) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        push @B, [map { -$_ } @$row];
    }

    __PACKAGE__->new(\@B);
}

sub map {
    my ($matrix, $callback) = @_;

    my @B;
    foreach my $row (@{$matrix->{A}}) {
        my @map;
        foreach my $elem (@$row) {
            local $_ = $elem;
            push @map, $callback->($elem);
        }
        push @B, \@map;
    }

    __PACKAGE__->new(\@B);
}

sub add {
    my ($m1, $m2) = @_;

    if (ref($m2) ne ref($m1)) {
        return $m1->scalar_add($m2);
    }

    my $A = $m1->{A};
    my $B = $m2->{A};

    my $rows = $#{$A};
    my $cols = $#{$A->[0]};

    my @C;
    foreach my $i (0 .. $rows) {
        push @C, [map { $A->[$i][$_] + $B->[$i][$_] } 0 .. $cols];
    }

    __PACKAGE__->new(\@C);
}

sub sub {
    my ($m1, $m2) = @_;

    my $r1 = ref($m1);
    my $r2 = ref($m2);

    if ($r1 ne $r2) {
        if ($r1 eq __PACKAGE__) {
            return $m1->scalar_sub($m2);
        }

        # a - b = -b + a
        if ($r2 eq __PACKAGE__) {
            return $m2->neg->scalar_add($m1);
        }
    }

    my $A = $m1->{A};
    my $B = $m2->{A};

    my $rows = $#{$A};
    my $cols = $#{$A->[0]};

    my @C;
    foreach my $i (0 .. $rows) {
        push @C, [map { $A->[$i][$_] - $B->[$i][$_] } 0 .. $cols];
    }

    __PACKAGE__->new(\@C);
}

sub mul {
    my ($m1, $m2) = @_;

    if (ref($m2) ne ref($m1)) {
        return $m1->scalar_mul($m2);
    }

    my $A = $m1->{A};
    my $B = $m2->{A};

    my @c;

    my $a_rows = $#{$A};
    my $b_rows = $#{$B};
    my $b_cols = $#{$B->[0]};

    foreach my $i (0 .. $a_rows) {
        foreach my $j (0 .. $b_cols) {
            foreach my $k (0 .. $b_rows) {

                my $t = $A->[$i][$k] * $B->[$k][$j];

                if (!defined($c[$i][$j])) {
                    $c[$i][$j] = $t;
                }
                else {
                    $c[$i][$j] += $t;
                }
            }
        }
    }

    __PACKAGE__->new(\@c);
}

sub div {
    my ($m1, $m2) = @_;

    if (ref($m2) ne ref($m1)) {
        return $m1->scalar_div($m2);
    }

    $m1->mul($m2->inv);
}

sub floor {
    my ($self) = @_;

    __PACKAGE__->new(
        [
         map {
             [
              map {
                  my $t = CORE::int($_);
                  $t -= 1 if ($_ != $t and $_ < 0);
                  $t;
                } @$_
             ]
         } @{$self->{A}}
        ]
    );
}

sub ceil {
    my ($self) = @_;

    __PACKAGE__->new(
        [
         map {
             [
              map {
                  my $t = CORE::int($_);
                  $t += 1 if ($_ != $t and $_ > 0);
                  $t;
                } @$_
             ]
         } @{$self->{A}}
        ]
    );
}

sub mod {
    my ($A, $B) = @_;

    my $r1 = ref($A);
    my $r2 = ref($B);

    if ($r1 ne $r2) {

        if ($r1 eq __PACKAGE__) {
            return $A->scalar_mod($B);
        }

        # A - B*floor(A/B)
        if ($r2 eq __PACKAGE__) {
            return Math::MatrixLUP::sub($A, $B->mul($B->map(sub { $A / $_ })->floor));
        }
    }

    # A - B*floor(A/B)
    $A->sub($B->mul($A->div($B)->floor));
}

sub pow {
    my ($A, $pow) = @_;

    $pow = CORE::int($pow);

    my $neg = 0;

    if ($pow < 0) {
        $neg = 1;
        $pow = -$pow;
    }

#<<<
        my $n = $#{$A->{A}};
        my $B = __PACKAGE__->new([map {
            my $i = $_;
            [map { ($i == $_) ? 1 : 0 } 0 .. $n]
        } 0 .. $n]);
#>>>

    return $B if ($pow == 0);

    while (1) {
        $B = $B->mul($A) if ($pow & 1);
        $pow >>= 1 or last;
        $A = $A->mul($A);
    }

    $neg ? $B->inv : $B;
}

sub solve {
    my ($self, $vector) = @_;

    my ($N, $A, $P) = @{$self->decompose};

    my @x = map { $vector->[$P->[$_]] } 0 .. $N;

    foreach my $i (1 .. $N) {
        foreach my $k (0 .. $i - 1) {
            $x[$i] -= $A->[$i][$k] * $x[$k];
        }
    }

    for (my $i = $N ; $i >= 0 ; --$i) {
        foreach my $k ($i + 1 .. $N) {
            $x[$i] -= $A->[$i][$k] * $x[$k];
        }
        $x[$i] /= $A->[$i][$i];
    }

    \@x;
}

sub invert {
    my ($self) = @_;

    my ($N, $A, $P) = @{$self->decompose};

    my @I;

    foreach my $j (0 .. $N) {
        foreach my $i (0 .. $N) {

            $I[$i][$j] = ($P->[$i] == $j) ? 1 : 0;

            foreach my $k (0 .. $i - 1) {
                $I[$i][$j] -= $A->[$i][$k] * $I[$k][$j];
            }
        }

        for (my $i = $N ; $i >= 0 ; --$i) {

            foreach my $k ($i + 1 .. $N) {
                $I[$i][$j] -= $A->[$i][$k] * $I[$k][$j];
            }

            $I[$i][$j] /= $A->[$i][$i] // return __PACKAGE__->new([[]]);
        }
    }

    __PACKAGE__->new(\@I);
}

*inv = \&invert;

sub determinant {
    my ($self) = @_;

    my ($N, $A, $P) = @{$self->decompose};

    my $det = $A->[0][0] // return 1;

    foreach my $i (1 .. $N) {
        $det *= $A->[$i][$i];
    }

    if (($P->[$N + 1] - $N) % 2 == 0) {
        $det *= -1;
    }

    $det;
}

*det = \&determinant;

sub stringify {
    my ($self) = @_;
    "[\n  " . join(",\n  ", map { "[" . join(", ", @$_) . "]" } @{$self->{A}}) . "\n]";
}

1;    # End of Math::MatrixLUP
