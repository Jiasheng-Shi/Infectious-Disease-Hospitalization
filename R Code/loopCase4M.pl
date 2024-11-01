## Odays-ellboot, Odays=114, 69 or 54
foreach my $i (1...10) {
  system "echo \"#!/bin/sh\" > settingM_$i\.sh";
  # system "echo \"BSUB -u Jiasheng.Shi\@pennmedicine.upenn.edu\" >> R0_covid_$i\.sh";
  system "echo \"module load R\" >> settingM_$i\.sh";

  # system "echo \"Rscript exec\.R $i \" >> setting_$i\.sh" ;
  system "echo \"Rscript simulation_hospitalization_organized_Case4\.R $i \" >> settingM_$i\.sh" ;
  system "bsub -e EM$i\.err -o OM$i\.out < settingM_$i\.sh";

}
