# LaTeX2HTML 96.1 (Feb 5, 1996)
# Associate images original text with physical files.

$key = q/{_inline}$calL${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=11 HEIGHT=13 ALIGN=BOTTOM ALT="tex2html_wrap_inline867" SRC="img12.gif"  > '; 
$key = q/{_inline}$kappa${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=8 HEIGHT=8 ALIGN=BOTTOM ALT="tex2html_wrap_inline943" SRC="img17.gif"  > '; 
$key = q/{_inline}$calP${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=12 HEIGHT=13 ALIGN=BOTTOM ALT="tex2html_wrap_inline863" SRC="img11.gif"  > '; 
$key = q/{_inline}$hats${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=7 HEIGHT=12 ALIGN=BOTTOM ALT="tex2html_wrap_inline785" SRC="img3.gif"  > '; 
$key = q/{_inline}$S_1,ldots,S_n${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=75 HEIGHT=25 ALIGN=MIDDLE ALT="tex2html_wrap_inline793" SRC="img4.gif"  > '; 
$key = q/{figure}framebox[5.75in]minipage5.25invskip.1inrmbfBranchingOperationunderbarInput:Asubproblem_inline$calS$_inlineand_inline$hatx$_inline,theLPsolutionyieldingthelowerbound.underbarOutput:_inline$S_1,ldots,S_p$_inlinesuchthat_inline$calS=cup_i=1^pS_i$_inline.bfStep1.Determinesets_inline$calL_1,ldots,calL_p$_inlineofinequalitiessuchthat_inline$calS=cup_i=1^n{xincalS:axleqbeta;forall(a,beta)incalL_i}$_inlineand_inline$hatxnotincup_i=1^nS_i$_inline.bfStep2.Set_inline$S_i={xincalS:axleqbeta;;forall;(a,beta)incalL_icupcalL'}$_inlinewhere_inline$calL'$_inlineisthesetofinequalitiesusedtodescribe_inline$calS$_inline.minipagelabelbranching-fig{figure}/;
$cached_env_img{$key} = ' <IMG WIDTH=661 HEIGHT=168 ALIGN=BOTTOM ALT="figure135" SRC="img14.gif"  > '; 
$key = q/{figure}framebox[5.75in]minipage5.25invskip.1inrmbfGenericBranchandCutAlgorithmunderbarInput:Adataarrayspecifyingtheprobleminstance.underbarOutput:Theglobaloptimalsolution_inline$s^*$_inlinetotheprobleminstance.bfStep1.Generatea``good''feasiblesolution_inline$hats$_inlineusingheuristics.Set_inline$alphaleftarrowc(hats)$_inline.bfStep2.Generatethefirstsubproblem_inline$calS^I$_inlinebyconstructingasmallset_inline$calL'$_inlineofinequalitiesvalidfor_inline$calP$_inline.Set_inline$Aleftarrow{calS^I}$_inline.bfStep3.If_inline$A=emptyset$_inline,STOPandoutput_inline$hats$_inlineastheglobaloptimum_inline$s^*$_inline.Otherwise,choosesome_inline$calSinA$_inline.Set_inline$AleftarrowAsetminus{calS}$_inline.Process_inline$calS$_inline.bfStep4.IftheresultofStep3isafeasiblesolution_inline$overlines$_inline,then_inline$coverlines;SPMlt;chats$_inline.Set_inline$hatsleftarrowoverlines$_inlineand_inline$alphaleftarrowc(overlines)$_inlineandgotoStep3.Ifthesubproblemwaspruned,gotoStep3.Otherwise,gotoStep5.bfStep5.Performthebranchingoperation.Addthesetofsubproblemsgeneratedto_inline$A$_inlineandgotoStep3.minipagelabelgb;SPMamp;c{figure}/;
$cached_env_img{$key} = ' <IMG WIDTH=661 HEIGHT=298 ALIGN=BOTTOM ALT="figure167" SRC="img15.gif"  > '; 
$key = q/{_inline}$cinbfR^S${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=51 HEIGHT=33 ALIGN=MIDDLE ALT="tex2html_wrap_inline775" SRC="img1.gif"  > '; 
$key = q/{_inline}$cup_i=1^nS_i=S${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=83 HEIGHT=24 ALIGN=MIDDLE ALT="tex2html_wrap_inline795" SRC="img5.gif"  > '; 
$key = q/{_inline}$calF${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=14 HEIGHT=13 ALIGN=BOTTOM ALT="tex2html_wrap_inline861" SRC="img10.gif"  > '; 
$key = q/{_inline}$calFsubseteq2^E${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=56 HEIGHT=33 ALIGN=MIDDLE ALT="tex2html_wrap_inline857" SRC="img8.gif"  > '; 
$key = q/{_inline}$cinbfR^E${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=52 HEIGHT=31 ALIGN=MIDDLE ALT="tex2html_wrap_inline859" SRC="img9.gif"  > '; 
$key = q/{_inline}$^copyright${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=15 HEIGHT=17 ALIGN=BOTTOM ALT="tex2html_wrap_inline945" SRC="img18.gif"  > '; 
$key = q/{figure}framebox[5.75in]minipage5.25invskip.1inrmbfBoundingOperationunderbarInput:Asubproblem_inline$calS$_inline,describedintermsofa``small''setofinequalities_inline$calL'$_inlinesuchthat_inline$calS={x^s:sincalF;hboxrmand;ax^sleqbeta;forall(a,beta)incalL'}$_inlineand_inline$alpha$_inline,anupperboundontheglobaloptimalvalue.underbarOutput:Either(1)anoptimalsolution_inline$s^*incalS$_inlinetothesubproblem,(2)alowerboundontheoptimalvalueofthesubproblem,or(3)amessagettprunedindicatingthatthesubproblemshouldnotbeconsideredfurther.bfStep1.Set_inline$calCleftarrowcalL'$_inline.bfStep2.SolvetheLP_inline$min{cx:axleqbeta;forall(a,beta)incalC}$_inline.bfStep3.IftheLPhasafeasiblesolution_inline$hatx$_inline,thengotoStep4.Otherwise,STOPandoutputttpruned.Thissubproblemhasnofeasiblesolutions.bfStep4.If_inline$chatx;SPMlt;alpha$_inline,thengotoStep5.Otherwise,STOPandoutputttpruned.Thissubproblemcannotproduceasolutionofvaluebetterthan_inline$alpha$_inline.bfStep5.If_inline$hatx$_inlineistheincidencevectorofsome_inline$hatsincalS$_inline,then_inline$hats$_inlineistheoptimalsolutiontothissubproblem.STOPandoutput_inline$hats$_inlineas_inline$s^*$_inline.Otherwise,applyseparationalgorithmsandheuristicsto_inline$hatx$_inlinetogetasetofviolatedinequalities_inline$calC'$_inline.If_inline$calC'=emptyset$_inline,then_inline$chatx$_inlineisalowerboundonthevalueofanoptimalelementof_inline$calS$_inline.STOPandreturn_inline$hatx$_inlineandthelowerbound_inline$chatx$_inline.Otherwise,set_inline$calCleftarrowcalCcupcalC'$_inlineandgotoStep2.minipagelabelproc-bound{figure}/;
$cached_env_img{$key} = ' <IMG WIDTH=661 HEIGHT=406 ALIGN=BOTTOM ALT="figure70" SRC="img6.gif"  > '; 
$key = q/{_inline}$hatsinS${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=39 HEIGHT=25 ALIGN=MIDDLE ALT="tex2html_wrap_inline777" SRC="img2.gif"  > '; 
$key = q/{_inline}$hboxemCP=(E,calF)${_inline}/;
$cached_env_img{$key} = ' <IMG WIDTH=94 HEIGHT=27 ALIGN=MIDDLE ALT="tex2html_wrap_inline853" SRC="img7.gif"  > '; 
$key = q/{figure}centeringpsfigfigure=hometkrPaperspicspbandc.pslabeloverview{figure}/;
$cached_env_img{$key} = ' <IMG WIDTH=633 HEIGHT=823 ALIGN=BOTTOM ALT="figure259" SRC="img16.gif"  > '; 

1;

