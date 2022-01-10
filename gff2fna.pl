#! /usr/bin/perl
############################################################
#       Copyright (C) 2021
#       Written by 葛文龙 (gwl9505@163.com)
#       Version: 1.1 (Sep 26th, 2021)
############################################################
use v5.12;

use Getopt::Long;
use FindBin qw($RealBin $Script);
use Cwd;
use List::MoreUtils ':all';
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my $line="————" x 20;
my $pl_site=cwd();

sub select{
        say "$line";
        say "\t-h  : 展示说明列表";
        say "\t  程序会按照gff给的位置信息抓出对应的fna序列\n";
        say "\t  使用举例，-g (gff注释文件) -f (基因组fna文件) -o (输出文件夹) -t";
				say "\t  特别说明，当启用-t参数时，会列出类型列表，选择输入所需要的类型即可\n";
        say "\t-g  : gff 格式的注释文件";
        say "\t-f  : 基因组数据文件\n";
        say "\t-o  : 输出文件夹名,默认为out_gff文件名";
				say "\t-t  : 类型筛选（不启用则默认输出全部）\n";
        say "$line";
        exit(1);
}
use vars qw($op_h $op_f $op_g $op_o $op_t);
GetOptions(
	"h" => \$op_h,
	"f:s" => \$op_f,
	"g:s" => \$op_g,
	"o:s" => \$op_o,
	"t" => \$op_t
);

&select if ($op_h);

##参数合法性检查
if(!$op_g){
	say $line;
	die "-g参数为空值，请检查\n";
}
if(!$op_f){
	say $line;
	die "-f参数为空值，请检查\n";
}
my $gff_site = $pl_site ."/$op_g";
my $fna_site = $pl_site ."/$op_f";

if(!-s $gff_site){
	die "指定的gff文件不存在或为空，请检查-g参数\n";
}
if(!-s $fna_site){
	die "指定的fna文件不存在或为空，请检查-f参数\n";
}
my $gff_name;
if($op_g =~/(.*)\./){
	$gff_name = $1;
}
$op_o ||= "gff2fa_out_$gff_name" ;
my $out_site = $pl_site ."/$op_o";

my $old_name= $out_site;

if(-e "$out_site"){
        my $num=1;
        $out_site.="_$num";
}

until(!-e "$out_site"){
        if($out_site=~/_([\d]+)$/){
        	my $n=$1 + 1;
      	  $out_site = $old_name."_$n";
        }
}

##主目录创建
say "\n开始运行，以下为提示信息";
say $line;
`mkdir $out_site`;

##gff类型解析
my @op_t;
my $type;
my @gff_type;
if($op_t){
	open GFF,"$gff_site";
	my @gff_type_list;
	while(<GFF>){
		chomp;
		next if(/\A#/);
		my @hang=split/\t/,$_;
		push @gff_type_list,$hang[2];
	}
	close GFF;
	@gff_type=uniq(@gff_type_list);
	say "\t类型列表如下，请选择需要输出的类型\n";
	say "\t输入 n ，则可以选择多个类型并放入对应文件夹中";
	say "\t输入 y ，则会将所有类型放进对应文件夹中：";
	say $line;
	foreach (@gff_type){say "\t$_ ";}
	print "\n请在此输入:\n";
	RE:chomp($op_t=<>);
	if($op_t eq "y"){
		$type="全";
		say BOLD GREEN "$type类型模式";
	}elsif($op_t eq "n"){
		$type="多";
		say BOLD GREEN "$type类型模式，输入想要的类型，按Ctrl+D结束";
		chomp (@op_t=<>);
	}elsif(grep /\b$op_t\b/,@gff_type){
		$type="单";
		say BOLD GREEN "$type类型模式";
	}else{
		say BOLD YELLOW "\n错误，输入的类型不在列表中或不是y/n，检查一下重新输入";
		goto RE;
	}
}

##fna处理
$/=">";
open FNA,$fna_site;
my %fna_gene;
while(<FNA>){
	chomp;
	my @s=split(/\n/,$_);
	my $n=$s[0];
	$n=$1 if($s[0]=~/^(.+?)(\s|\t)/);
	if(defined($n)){
		my $s=join("",@s[1..$#s]);
		$fna_gene{$n}=$s;
	}else {
		say "$n";
	}
}
close $fna_site;
$/="\n";

##gff信息提取
my $star_num;
my $end_num;
my $id;

my @gff_gene_list_all;

open GFF,"$gff_site";
while(<GFF>){
	next if(/\A#/);
	my @hang=split/\t/,$_;
	push @gff_gene_list_all,$hang[0];
}
close GFF;
my @gff_gene_list=uniq(@gff_gene_list_all);
my @fna_gene_list=keys %fna_gene;
my $fna_gene_site;
my %id_type;
my $type_site;
for(@fna_gene_list){
		my %id_check;
		my $fna_check=$_;
    my $pip=grep/$fna_check/,@gff_gene_list;
		if($pip==1){
			$fna_gene_site=$out_site."/$fna_check";
			`mkdir $fna_gene_site`;
			open GFF,"$gff_site";
			my $old_fna_gene_site=$fna_gene_site;
			OP:while(<GFF>){
				chomp;
				next if(/\A#/);
				if(/\A$fna_check/){
					my @hang_cut=split/\t/,$_;
					$star_num = $hang_cut[4]-$hang_cut[3];
					my $xulie_len= length $fna_gene{$hang_cut[0]};
					if($star_num < 0){
						die "怪事，序列结束位置小于起始位置，我应该修正这个脚本以适应这么麻烦的gff\n程序中止\n";
					}
					my $hang=$_;
					if($type eq "单"){
						if($op_t ne $hang_cut[2]){goto OP;}
					}elsif($type eq "多"){
						my $multi_check=grep /\b$hang_cut[2]\b/,@op_t;
						if($multi_check == 0){goto OP;}
					}
					if($hang=~/ID=(.*?);/){
	        	       			$id=$1;
						$id_check{$id}+=1;
	          if($id_check{$id}>=2){
	          	print BOLD YELLOW "\n注意";
	           	print "，发现重复ID $id";
							$id.="_1";
							say "   已重命名为 $id\n";
	          }
						if($hang_cut[4] > $xulie_len){
							print BOLD RED "\n警告：";
							say " $id 所在序列区域超出 $fna_check 总序列长度，该ID已被跳过";
							goto REGE;
						}
						my $id_site;
						if($type eq "单"){
							$type_site=$fna_gene_site."/$op_t";
							if(!-s $type_site){`mkdir $type_site`;}
							$id_site=$type_site."/$id";
						}elsif ($type eq "多") {
							$id_type{$id}=$hang_cut[2];
							$type_site=$fna_gene_site."/$id_type{$id}";
							if(!-s $type_site){`mkdir $type_site`;}
							$id_site=$type_site."/$id";
						}elsif($type eq "全"){
							$type_site=$fna_gene_site."/$hang_cut[2]";
							if(!-s $type_site){`mkdir $type_site`;}
							$id_site=$type_site."/$id";
					}else{
						$id_site=$fna_gene_site."/$id";
					}
					open IN,">$id_site";
					my $num=$hang_cut[3]-1;
					my $xulie=substr($fna_gene{$hang_cut[0]},$num,$star_num);
					if($hang_cut[6] eq "-"){
						$xulie=~ tr/atcgATCG/tagcTAGC/;
	        				$xulie=reverse $xulie;
					}
					print IN "ID=$id\t$hang_cut[6]\t$hang_cut[3]\t$hang_cut[4]\n";
					print IN "$xulie";
					close IN;
					REGE:undef $star_num;
					undef $num;
					}
				}
			}
		}else{
			print BOLD YELLOW "\n注意";
			say "，$op_f 文件中读取到基因名 $_ ，在 $op_g 文件中不存在，已略过";
		}
			undef %id_check;
}
close GFF;
say "\n$line";
if($out_site=~/(.*)\//){
	say "已将结果保存至 $'";
}
say "程序结束，感谢使用~(*^_^*)";
