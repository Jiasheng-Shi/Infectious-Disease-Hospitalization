## Odays-ellboot, Odays=114, 69 or 54
foreach my $i (1...20) {
  system "echo \"#!/bin/sh\" > settingUO_$i\.sh";
  # system "echo \"BSUB -u Jiasheng.Shi\@pennmedicine.upenn.edu\" >> R0_covid_$i\.sh";
  system "echo \"module load R\" >> settingUO_$i\.sh";

  # system "echo \"Rscript exec\.R $i \" >> setting_$i\.sh" ;
  system "echo \"Rscript simulation_hospitalization_organized_Case3UO\.R $i \" >> settingUO_$i\.sh" ;
  system "bsub -e EUO$i\.err -o OUO$i\.out < settingUO_$i\.sh";

}
