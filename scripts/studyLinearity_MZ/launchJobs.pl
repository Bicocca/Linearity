use Env;



open (USERCONFIG,$ARGV[0]) ;
while (<USERCONFIG>)
{
  chomp; 
  s/#.*//;                # no comments
  s/^\s+//;               # no leading white
  s/\s+$//;               # no trailing white
  my ($var, $value) = split(/\s*=\s*/, $_, 2);
  $User_Preferences{$var} = $value;
}

$JOBCfgTemplate      = $User_Preferences{"JOBCfgTemplate"};
$INPUTFilesDA        = $User_Preferences{"INPUTFilesDA"};
$INPUTFilesMC        = $User_Preferences{"INPUTFilesMC"};
$USEGlobeNtuple      = $User_Preferences{"USEGlobeNtuple"};
$USEShervinNtuple    = $User_Preferences{"USEShervinNtuple"};
$MCGen               = $User_Preferences{"MCGen"};
$RUNDepFlag          = $User_Preferences{"RUNDepFlag"};
$ENCorrType          = $User_Preferences{"ENCorrType"};
$YEAR                = $User_Preferences{"YEAR"};
$DATALabel           = $User_Preferences{"DATALabel"};
$CATType             = $User_Preferences{"CATType"};
$EVTSPerPoint        = $User_Preferences{"EVTSPerPoint"};
$ENERGYScaleCorrType = $User_Preferences{"ENERGYScaleCorrType"};
$ENERGYSmearingType  = $User_Preferences{"ENERGYSmearingType"};
$MCClosure           = $User_Preferences{"MCClosure"};
$MCHiggs             = $User_Preferences{"MCHiggs"};
$DPHIMax             = $User_Preferences{"DPHIMax"};
$NCats               = $User_Preferences{"NCats"};



$baseDir = $ENV{LINEARITY};
$DphiLabel = sprintf("Dphi%dp%02d",int($DPHIMax),int($DPHIMax*100)%100);


$globeLabel = "globe";
if( $USEGlobeNtuple eq "false" )
{
  $globeLabel = "nonGlobe";
}

$runDepLabel = "allRange";
if( $RUNDepFlag eq "true" )
{
  $runDepLabel = "runDependent";
}

$mcClosureLabel = "";
if( $MCClosure eq "true" )
{
  $mcClosureLabel = "_MCClosure";
}

$mcHiggsLabel = "";
if( $MCHiggs eq "true" )
{
  $mcHiggsLabel = "_MCHiggs";
}

$jobDir = $baseDir."/data/".$YEAR."/".$DATALabel."/MZ_".$CATType."_".$globeLabel."_".$MCGen."-".$runDepLabel."_".$ENCorrType."_".$DphiLabel.$mcClosureLabel.$mcHiggsLabel."/jobs/";
print("jobDir: ".$jobDir."\n");

$command = "mkdir -p ".$jobDir;
system($command);


@evtsPerPoint = split(',',$EVTSPerPoint);


open(LANCIA,">","lancia.sh") or die "Can't open file lancia.sh";


for($cat = 0; $cat < $NCats; ++$cat)
{
  print(">>> cat: ".$cat."\n");
  
  
  $jobCfg = $jobDir."/params_cat".$cat.".cfg";
  system ("cat ".$JOBCfgTemplate."   | sed -e s%INPUTFILESDA%".$INPUTFilesDA.
                                 "%g | sed -e s%INPUTFILESMC%".$INPUTFilesMC.
                                 "%g | sed -e s%USEGLOBENTUPLE%".$USEGlobeNtuple.
                                 "%g | sed -e s%MCGEN%".$MCGen.
                                 "%g | sed -e s%RUNDEPFLAG%".$RUNDepFlag.
                                 "%g | sed -e s%USESHERVINNTUPLE%".$USEShervinNtuple.
                                 "%g | sed -e s%ENCORRTYPE%".$ENCorrType.
                                 "%g | sed -e s%YEAR%".$YEAR.
                                 "%g | sed -e s%DATALABEL%".$DATALabel.
                                 "%g | sed -e s%CATTYPE%".$CATType.
                                 "%g | sed -e s%CATEGORY%".$cat.
                                 "%g | sed -e s%EVTSPERPOINT%".$evtsPerPoint[$cat].
                                 "%g | sed -e s%ENERGYSCALECORRTYPE%".$ENERGYScaleCorrType.
                                 "%g | sed -e s%ENERGYSMEARINGTYPE%".$ENERGYSmearingType.
                                 "%g | sed -e s%MCCLOSURE%".$MCClosure.
                                 "%g | sed -e s%MCHIGGS%".$MCHiggs.
                                 "%g | sed -e s%DPHIMAX%".$DPHIMax.
                                 "%g > ".$jobCfg);
  
  
  $jobOut = $jobDir."/out_cat".$cat.".txt";
  
  $jobSh = $jobDir."/job_cat".$cat.".sh";
  open(JOBSH,">",$jobSh) or die "Can't open file ".$jobSh;
  
  print JOBSH "echo \$SHELL\n";
  print JOBSH "source /gwpool/initcms/slc5_64-gcc462-cmssw.sh\n";
  print JOBSH "cd ".$baseDir."\n";
  print JOBSH "cd ../CMSSW_5_3_5/src\n";
  print JOBSH "eval `scramv1 runtime -sh`\n";
  print JOBSH "cd -\n";
  print JOBSH "pwd \n";
  print JOBSH "source scripts/setup.sh \n";
  print JOBSH "root-config --version \n";
  print JOBSH "echo \$PATH\n";
  print JOBSH "unbuffer studyLinearity_MZ.exe ".$jobCfg." >& ".$jobOut."\n";
  
  
  print LANCIA "qsub -d ".$jobDir." -q shortcms ".$jobSh."\n";
}
