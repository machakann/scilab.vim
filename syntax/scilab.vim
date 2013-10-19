" Vim syntax file
" Language:	Scilab
" Maintainer:	mckn <mckn@outlook.jp> from Patricio Toledo <patoledo@ing.puc.cl>, from matlab.sci
" Last Change: 19-Oct-2013.

" For version 5.x: Clear all syntax items
" For version 6.x: Quit when a syntax file was already loaded
if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif

" Scilab uses % as a part of function name and constant name.
setl iskeyword+=%

syn keyword scilabStatement		function endfunction return abort pause resume
syn keyword scilabLabel			case switch select
syn keyword scilabConditional		else elseif end if otherwise then
syn keyword scilabRepeat		do for while break continue
syn keyword scilabOperator		and or
syn keyword scilabSpecialVariables	varargin varargout
syn keyword scilabPreDefVariables	SCI SCIHOME TMPDIR %eps %i %inf %nan %pi
syn keyword scilabException		try catch

syn keyword scilabTodo			contained  TODO

" If you do not want these operators lit, uncommment them and the "hi link" below
syn match scilabArithmeticOperator	"[-+]"
syn match scilabArithmeticOperator	"\.\=[*/\\^]"
syn match scilabRelationalOperator	"[=~]="
syn match scilabRelationalOperator	"[<>]=\="
syn match scilabLogicalOperator		"[&|~]"
syn match scilabLineContinuation	"\.\{2}\ze.*$"

"syn match scilabIdentifier		"\<\a\w*\>"

" String
syn region scilabString			start=+'+ end=+'+	oneline

" If you don't like tabs
syn match scilabTab			"\t"
" Standard numbers
syn match scilabNumber			"\<\d\+[ij]\=\>"
" floating point number, with dot, optional exponent
syn match scilabFloat			"\<\d\+\(\.\d*\)\=\([edED][-+]\=\d\+\)\=[ij]\=\>"
" floating point number, starting with a dot, optional exponent
syn match scilabFloat			"\.\d\+\([edED][-+]\=\d\+\)\=[ij]\=\>"
" Transpose character and delimiters: Either use just [...] or (...) aswell
syn match scilabDelimiter		"[][]"
syn match scilabDelimiter		"[][()]"
syn match scilabTransposeOperator	"[])a-zA-Z0-9.]'"lc=1
syn match scilabSizeOperator		"\$"
syn match scilabSemicolon		";"
syn match scilabColon			","
syn match scilabComment			"//.*$"	contains=scilabTodo,scilabTab
syn match scilabError			"-\=\<\d\+\.\d\+\.[^*/\\^]"
syn match scilabError			"-\=\<\d\+\.\d\+[eEdD][-+]\=\d\+\.\([^*/\\^]\)"

syn keyword scilabScilab		banner getdebuginfo getmemory getmodules
syn keyword scilabScilab		getos getscilabmode getshell getversion
syn keyword scilabScilab		gstacksize stacksize ver with_javasci
syn keyword scilabScilab		with_macros_source with_module with_tk
syn keyword scilabScilab		debug delbpt dispbpt setbpt where whereami
syn keyword scilabScilab		errcatch errclear error error_table
syn keyword scilabScilab		iserror lasterror warning
syn keyword scilabScilab		clear clearglobal exists getvariablesonstack
syn keyword scilabScilab		isdef isglobal names predef who who_user whos
syn keyword scilabScilab		exit perl quit scilab startup
syn keyword scilabDiffAndIntegral	bvodeS dae daeoptions dasrt dassl diff
syn keyword scilabDiffAndIntegral	feval impl int2d int3d intc integrate
syn keyword scilabDiffAndIntegral	intg intl intsplin inttrap numdiff ode
syn keyword scilabDiffAndIntegral	ode_discrete ode_optional_output
syn keyword scilabDiffAndIntegral	ode_root odedc odeoptions
syn keyword scilabElementaryFunction	bitand bitcmp bitget bitor
syn keyword scilabElementaryFunction	bitset bitxor isequalbitwise
syn keyword scilabElementaryFunction	complex conj imag imult isreal real
syn keyword scilabElementaryFunction	factor factorial gcd lcm perms primes rat
syn keyword scilabElementaryFunction	diag eye ind2sub linspace logspace
syn keyword scilabElementaryFunction	meshgrid ndgrid ones rand squarewave
syn keyword scilabElementaryFunction	sub2ind testmatrix toeplitz zeros
syn keyword scilabElementaryFunction	exp expm log log10 log1p log2
syn keyword scilabElementaryFunction	logm polar sqrt sqrtm
syn keyword scilabElementaryFunction	ceil clean double fix floor format
syn keyword scilabElementaryFunction	frexp ieee int isinf isnan nearfloat
syn keyword scilabElementaryFunction	nextpow2 number_properties round
syn keyword scilabElementaryFunction	base2dec bin2dec dec2base dec2bin
syn keyword scilabElementaryFunction	dec2hex dec2oct hex2dec oct2dec
syn keyword scilabElementaryFunction	flipdim matrix permute pertrans
syn keyword scilabElementaryFunction	repmat resize_matrix squeeze
syn keyword scilabElementaryFunction	abs cumprod cumsum kron max min
syn keyword scilabElementaryFunction	norm prod signm sum tril triu
syn keyword scilabElementaryFunction	dsearch gsort lex_sort vectorfind
syn keyword scilabElementaryFunction	intersect setdiff union unique
syn keyword scilabElementaryFunction	bloc2exp bloc2ss pen2ea ssrand
syn keyword scilabElementaryFunction	sysconv sysdiag trfmod
syn keyword scilabElementaryFunction	addf cmb_lin ldivf mulf rdivf
syn keyword scilabElementaryFunction	solve subf trianfml trisolve
syn keyword scilabElementaryFunction	acos acosd acosh acoshm acosm acot
syn keyword scilabElementaryFunction	acotd acoth acsc acscd acsch asec
syn keyword scilabElementaryFunction	asecd asech asin asind asinh asinhm
syn keyword scilabElementaryFunction	asinm atan atand atanh atanhm atanm
syn keyword scilabElementaryFunction	cos cosd cosh coshm cosm cotd cotg
syn keyword scilabElementaryFunction	coth cothm csc cscd csch csgn sec
syn keyword scilabElementaryFunction	secd sech sin sinc sind sinh sinhm
syn keyword scilabElementaryFunction	sinm tan tand tanh tanhm tanm
syn keyword scilabElementaryFunction	cat cell2mat cellstr isempty
syn keyword scilabElementaryFunction	isequal isvector lstsize pmodulo
syn keyword scilabElementaryFunction	ndims nthroot sign size
syn keyword scilabLinearAlgebra		balanc bdiag gschur gspec hess pbig
syn keyword scilabLinearAlgebra		projspec psmall schur spec sva svd
syn keyword scilabLinearAlgebra		givens householder sqroot
syn keyword scilabLinearAlgebra		colcomp fullrf fullrfk im_inv
syn keyword scilabLinearAlgebra		kernel range rowcomp
syn keyword scilabLinearAlgebra		aff2ab chol inv linsolve
syn keyword scilabLinearAlgebra		lsq lu pinv qr rankqr
syn keyword scilabLinearAlgebra		classmarkov eigenmarkov genmarkov
syn keyword scilabLinearAlgebra		cond det orth rank rcond rref trace
syn keyword scilabLinearAlgebra		companion ereduc fstair glever kroneck
syn keyword scilabLinearAlgebra		lyap pencan penlaur quaskro randpencil
syn keyword scilabLinearAlgebra		rowshuff sylv
syn keyword scilabLinearAlgebra		coff nlev
syn keyword scilabLinearAlgebra		spaninter spanplus spantwo
syn keyword scilabLinearAlgebra		proj
syn keyword scilabInterpolation		bsplin3val cshep2d eval_cshep2d interp interp1
syn keyword scilabInterpolation		interp2d interp3d interpln linear_interpn lsq_splin
syn keyword scilabInterpolation		smooth splin splin2d splin3d
syn keyword scilabCACSD			abcd cont_frm dbphi des2ss des2tf frep2tf lsslist
syn keyword scilabCACSD			markp2ss sm2des sm2ss ss2des ss2ss ss2tf tf2des tf2ss
syn keyword scilabCACSD			black bode chart evans gainplot hallchart m_circle
syn keyword scilabCACSD			nicholschart nyquist phaseplot sgrid show_margins svplot zgrid
syn keyword scilabCACSD			abinv arhnk arl2 arma arma2p arma2ss armac armax armax1
syn keyword scilabCACSD			arsimul augment balreal bilin bstap cainv calfrq canon
syn keyword scilabCACSD			ccontrg cls2dls colinout colregul cont_mat contr contrss
syn keyword scilabCACSD			copfac csim ctr_gram damp dcf ddp dhinf dhnorm dscr dsimul
syn keyword scilabCACSD			dt_ility dtsi equil equil1 feedback findABCD findAC findBD
syn keyword scilabCACSD			findBDK findR findx0BD flts fourplan freq freson fspec fspecg
syn keyword scilabCACSD			fstabst g_margin gamitg gcare gfare gfrancis gtild h2norm h_cl
syn keyword scilabCACSD			h_inf h_inf_st h_norm hankelsv hinf imrep2ss inistate invsyslin
syn keyword scilabCACSD			kpure krac2 lcf leqr lft lin linf linfn linmeq lqe lqg lqg2stan
syn keyword scilabCACSD			lqg_ltr lqr ltitr macglov minreal minss mucomp narsimul nehari
syn keyword scilabCACSD			noisegen nyquistfrequencybounds obs_gram obscont observer obsv_mat
syn keyword scilabCACSD			obsvss p_margin parrot pfss phasemag plzr pol2des ppol prbs_a projsl
syn keyword scilabCACSD			reglin repfreq ric_desc ricc riccati routh_t rowinout rowregul rtitr
syn keyword scilabCACSD			sensi sident sorder specfact ssprint st_ility stabil sysfact syslin
syn keyword scilabCACSD			syssize time_id trzeros ui_observer unobs zeropen
syn keyword scilabPolynomials		bezout chepol cmndred coeff coffg colcompr degree denom derivat
syn keyword scilabPolynomials		determ detr diophant factors hermit horner hrmt htrianr inv_coeff
syn keyword scilabPolynomials		invr lcmdiag ldiv numer pdiv pol2str polfact poly rational residu
syn keyword scilabPolynomials		roots rowcompr sfact simp simp_mode sylm systmat varn
syn keyword scilabSignalProcessing	analpf buttmag casc cheb1mag cheb2mag convol ell1mag eqfir eqiir
syn keyword scilabSignalProcessing	faurre ffilt filter find_freq frmag fsfirlin group iir iirgroup
syn keyword scilabSignalProcessing	iirlp kalm lev levin lindquist remezb srfaur srkf sskf system trans
syn keyword scilabSignalProcessing	wfir wiener wigner window zpbutt zpch1 zpch2 zpell
syn keyword scilabSignalProcessing	lattn lattp phc rpem
syn keyword scilabSignalProcessing	bilt jmat sincd
syn keyword scilabSignalProcessing	corr cspect czt intdec mese pspect
syn keyword scilabSignalProcessing	idct dft (deprecated) idst ifft hank hilb mfft
syn keyword scilabSignalProcessing	cepstrum conv conv2 convol2d detrend fft2 fftshift filt_sinc frfit
syn keyword scilabSignalProcessing	hilbert mrfit remez syredi wfir_gui xcorr xcov yulewalk
syn keyword scilabFFTW			fftw_flags fftw_forget_wisdom get_fftw_wisdom set_fftw_wisdom
syn keyword scilabSpecialFunctions	amell besselh beta calerf delip dlgamma erf erfc erfcx
syn keyword scilabSpecialFunctions	erfinv findm gamma gammaln legendre %asn %k %sn
syn keyword scilabRandlib		grand
syn keyword scilabARnoldiPACKage	dnaupd dneupd dsaupd dseupd eigs znaupd zneupd
syn keyword scilabStatistics		binomial cdfbet cdfbin cdfchi cdfchn cdff cdffnc cdfgam cdfnbn cdfnor cdfpoi cdft
syn keyword scilabStatistics		geomean harmean mean meanf trimmean
syn keyword scilabStatistics		nancumsum nand2mean nanmax nanmean nanmeanf nanmedian nanmin nanstdev nansum thrownan
syn keyword scilabStatistics		center correl covar median msd mvvacov stdev stdevf variance variancef wcenter
syn keyword scilabStatistics		ftest ftuneq
syn keyword scilabStatistics		iqr mad strange
syn keyword scilabStatistics		cmoment moment perctl quart
syn keyword scilabStatistics		pca princomp show_pca
syn keyword scilabStatistics		regress
syn keyword scilabStatistics		sample samplef samwr
syn keyword scilabStatistics		nfreq tabul
syn keyword scilabSparseMatrix		ludel lufact luget lusolve spchol
syn keyword scilabSparseMatrix		gmres pcg qmr
syn keyword scilabSparseMatrix		issparse nnz speye spones sprand spzeros
syn keyword scilabSparseMatrix		adj2sp full mtlb_sparse sp2adj sparse spcompack spget
syn keyword scilabSparseMatrix		chfact chsolve ordmmd
syn keyword scilabUMFPACKInterface	PlotSparse ReadHBSparse cond2sp condestsp
syn keyword scilabUMFPACKInterface	rafiter res_with_prec taucs_chdel taucs_chfact
syn keyword scilabUMFPACKInterface	taucs_chget taucs_chinfo taucs_chsolve
syn keyword scilabUMFPACKInterface	umf_ludel umf_lufact umf_luget umf_luinfo
syn keyword scilabUMFPACKInterface	umf_lusolve umfpack
syn keyword scilabOptAndSim		fminsearch neldermead overview nmplot optimget
syn keyword scilabOptAndSim		optimplotfunccount optimplotfval optimplotx optimset
syn keyword scilabOptAndSim		datafit fit_dat leastsq lsqrsolve
syn keyword scilabOptAndSim		optimbase
syn keyword scilabOptAndSim		optimsimplex
syn keyword scilabOptAndSim		aplat list2vec lmisolver lmitool pack recons semidef unpack vec2list
syn keyword scilabOptAndSim		NDcost derivative fsolve karmarkar optim qld qp_solve qpsolve readmps
syn keyword scilabGeneticAlgorithms	optim_ga optim_moga optim_nsga optim_nsga2
syn keyword scilabGeneticAlgorithms	coding_ga_binary coding_ga_identity crossover_ga_binary crossover_ga_default
syn keyword scilabGeneticAlgorithms	init_ga_default mutation_ga_binary mutation_ga_default pareto_filter
syn keyword scilabGeneticAlgorithms	selection_ga_elitist selection_ga_random
syn keyword scilabSimulatedAnnealing	optim_sa
syn keyword scilabSimulatedAnnealing	accept_func_default accept_func_vfsa compute_initial_temp neigh_func_csa
syn keyword scilabSimulatedAnnealing	neigh_func_default neigh_func_fsa neigh_func_vfsa temp_law_csa
syn keyword scilabSimulatedAnnealing	temp_law_default temp_law_fsa temp_law_huang temp_law_vfsa
syn keyword scilabXMLManagement		XML Objects xmlAddNs xmlAppend xmlAsNumber xmlAsText xmlDTD xmlDelete
syn keyword scilabXMLManagement		xmlDocument xmlDump xmlElement xmlFormat xmlGetNsByHref xmlGetNsByPrefix
syn keyword scilabXMLManagement		xmlGetOpenDocs xmlIsValidObject xmlName xmlNs xmlRead xmlReadStr xmlRelaxNG
syn keyword scilabXMLManagement		xmlRemove xmlSchema xmlSetAttributes xmlValidate xmlWrite xmlXPath
syn keyword scilabFilesIO		cd createdir dir isdir ls mkdir pwd removedir rmdir
syn keyword scilabFilesIO		basename dirname fileext fileparts filesep fullfile fullpath
syn keyword scilabFilesIO		get_absolute_file_path getdrives getlongpathname getrelativefilename
syn keyword scilabFilesIO		getshortpathname is_absolute_path pathconvert pathsep tempname
syn keyword scilabFilesIO		copyfile deletefile dispfiles fileinfo findfiles fprintf fprintfMat fscanf
syn keyword scilabFilesIO		fscanfMat getmd5 isfile listfiles listvarinfile maxfiles mclearerr %io
syn keyword scilabFilesIO		mclose mdelete meof merror mfprintf msscanf mgeti mgetl mgetstr mopen movefile
syn keyword scilabFilesIO		mput mputl mputstr mseek mtell newest save format scanf scanf_conversion sscanf
syn keyword scilabIOFunctions		file getenv getio getpid getscilabkeywords halt host input load read read4b
syn keyword scilabIOFunctions		readb save setenv unix unix_g unix_s unix_w unix_x writb write write4b
syn keyword scilabGraphics		LineSpec Matplot Matplot1 Matplot properties Sfgrayplot Sgrayplot champ champ1
syn keyword scilabGraphics		champ properties comet contour2d contour2di contourf errbar fchamp fcontour2d
syn keyword scilabGraphics		fec fec properties fgrayplot fplot2d grayplot grayplot properties graypolarplot
syn keyword scilabGraphics		histplot paramfplot2d plot plot2d plot2d1 plot2d2 plot2d3 plot2d4 polarplot
syn keyword scilabGraphics		comet3d contour eval3d eval3dp fac3d fcontour fplot3d fplot3d1 genfac3d geom3d
syn keyword scilabGraphics		hist3d mesh milk_drop nf3d param3d param3d1 param3d properties plot3d plot3d1
syn keyword scilabGraphics		plot3d2 plot3d3 secto3d surf surface properties
syn keyword scilabGraphics		captions legend legends title xlabel ylabel zlabel xtitle
syn keyword scilabGraphics		gca gda graduate isoview newaxes plotframe replot rotate_axes
syn keyword scilabGraphics		sca sda square subplot unzoom zoom_rect
syn keyword scilabGraphics		drawaxis
syn keyword scilabGraphics		bar barh barhomogenize
syn keyword scilabGraphics		addcolor autumncolormap bonecolormap color color_list colorbar colordef
syn keyword scilabGraphics		colormap coolcolormap coppercolormap getcolor graycolormap hotcolormap
syn keyword scilabGraphics		hsv2rgb hsvcolormap jetcolormap name2rgb oceancolormap pinkcolormap
syn keyword scilabGraphics		rainbowcolormap rgb2name springcolormap summercolormap whitecolormap wintercolormap
syn keyword scilabGraphics		datatipCreate datatipGetEntities datatipGetStruct datatipInitStruct
syn keyword scilabGraphics		datatipManagerMode datatipMove datatipRedraw datatipRemove datatipRemoveAll
syn keyword scilabGraphics		datatipSetDisplay datatipSetInterp datatipSetOrientation datatipSetStruct
syn keyword scilabGraphics		datatipSetStyle datatipToggle datatips orthProj
syn keyword scilabGraphics		clf drawlater drawnow figure properties gcf gdf scf sdf
syn keyword scilabGraphics		xarc xarcs xarrows xfarc xfarcs xfrect xrect xrects
syn keyword scilabGraphics		copy delete draw gce ged get_figure_handle glue is_handle_valid relocate_handle swap_handles unglue
syn keyword scilabGraphics		dragrect edit_curv event_handler_functions locate rubberbox seteventhandler xclick xgetmouse
syn keyword scilabGraphics		xload xsave
syn keyword scilabGraphics		pie
syn keyword scilabGraphics		xfpoly xfpolys xpoly xpolys xrpoly
syn keyword scilabGraphics		get set
syn keyword scilabGraphics		getlinestyle getmark getsymbol
syn keyword scilabGraphics		getfont graphics stringbox titlepage xinfo xlfont xstring xstringb xstringl
syn keyword scilabGraphics		move rotate scaling
syn keyword scilabGraphics		havewindow show_window winsid
syn keyword scilabGraphics		clear_pixmap pixel_drawing_mode show_pixmap twinkle xchange xclear xdel
syn keyword scilabGraphics		xget xgetech xgraduate xgrid xname xnumb xpause xsegs xset xsetech xsetm
syn keyword scilabGraphics		driver xend xinit xs2bmp xs2emf xs2eps xs2gif xs2jpg xs2pdf xs2png xs2ppm xs2ps xs2svg
syn keyword scilabGUI			uiConcatTree uiCreateNode uiCreateTree uiDeleteNode uiDisplayTree uiDumpTree uiEqualsTree
syn keyword scilabGUI			uiFindNode uiGetChildrenNode uiGetNodePosition uiGetParentNode uiInsertNode
syn keyword scilabGUI			about addmenu clipboard close delmenu exportUI figure findobj gcbo getcallbackobject
syn keyword scilabGUI			getinstalledlookandfeels getlookandfeel getvalue messagebox printfigure printsetupbox
syn keyword scilabGUI			progressionbar root_properties setlookandfeel setmenu toolbar toprint tree_show uicontextmenu
syn keyword scilabGUI			uicontrol uigetcolor uigetdir uigetfile uigetfont uimenu uiputfile unsetmenu usecanvas waitbar
syn keyword scilabGUI			x_choices x_choose x_choose_modeless x_dialog x_matrix x_mdialog
syn keyword scilabDataStructure		cell definedfields fieldnames getfield hypermat hypermatrices iscell iscellstr isfield
syn keyword scilabDataStructure		isstruct list lstcat matrices mlist null rlist setfield struct tlist type typename typeof
syn keyword scilabParameters		add_param get_param init_param is_param list_param remove_param set_param
syn keyword scilabBoolean		bool2s find
syn keyword scilabIntegers		iconvert uint32 inttype
syn keyword scilabStrings		ascii asciimat blanks char convstr emptystr eval evstr grep isalphanum isascii isdigit isletter
syn keyword scilabStrings		isnum justify length part regexp sci2exp strcat strchr strcmp strcmpi strcspn strindex string
syn keyword scilabStrings		strings stripblanks strncpy strrchr strrev strsplit strspn strstr strsubst strtod strtok tokenpos tokens
syn keyword scilabSoundFileHandling	analyze auread auwrite beep lin2mu loadwave mapsound mu2lin playsnd savewave sound soundsec wavread wavwrite
syn keyword scilabTimeAndDate		calendar clock date datenum datevec eomday etime getdate now realtime sleep tic timer toc weekday
syn keyword scilabOutputFunctions	disp printf sprintf prettyprint print printf_conversion
" syn keyword scilabXcos			lincos scicos_simulate scicosim steadycos xcosValidateBlockSet xcosValidateCompareBlock xcos_simulate
" syn keyword scilabXcos			Annotations_pal TEXT_f
" syn keyword scilabXcos			Commonlyusedblocks_pal LOGICAL_OP RELATIONALOP
" syn keyword scilabXcos			Continuous_pal CLINDUMMY_f CLR CLSS DERIV INTEGRAL_f INTEGRAL_m PID TCLSS TIME_DELAY VARIABLE_DELAY
" syn keyword scilabXcos			Demonstrationsblocks_pal AUTOMAT BOUNCE BOUNCEXY BPLATFORM PDE
" syn keyword scilabXcos			discontinuities_pal BACKLASH DEADBAND HYSTHERESIS RATELIMITER SATURATION
" syn keyword scilabXcos			Discrete_pal DELAYV_f DELAY_f DLR DLRADAPT_f DLSS DOLLAR_f REGISTER
" syn keyword scilabXcos			Electrical_pal CCS CVS Capacitor ConstantVoltage CurrentSensor Diode Ground Gyrator IdealTransformer Inductor
" syn keyword scilabXcos			NMOS NPN OpAmp PMOS PNP PotentialSensor Resistor SineVoltage Switch VVsourceAC VariableResistor VoltageSensor VsourceAC
" syn keyword scilabXcos			Events_pal ANDBLK ANDLOG_f CEVENTSCOPE CLKFROM CLKGOTO CLKGotoTagVisibility CLKSOMV_f EDGE_TRIGGER ESELECT_f EVTDLY_c
" syn keyword scilabXcos			EVTGEN_f EVTVARDLY Extract_Activation HALT_f IFTHEL_f MCLOCK_f MFCLCK_f M_freq VirtualCLK0 freq_div
" syn keyword scilabXcos			Implicit_pal CONSTRAINT_c DIFF_f
" syn keyword scilabXcos			Integer_pal BITCLEAR BITSET CONVERT DFLIPFLOP DLATCH EXTRACTBITS INTMUL JKFLIPFLOP LOGIC SHIFT SRFLIPFLOP
" syn keyword scilabXcos			Lookuptables_pal INTRP2BLK_f INTRPLBLK_f LOOKUP_f
" syn keyword scilabXcos			Mathoperations_pal ABS_VALUE BIGSOM_f COSBLK_f EXPBLK_m GAINBLK_f INVBLK LOGBLK_f MATMAGPHI MATZREIM
" syn keyword scilabXcos			MAXMIN MAX_f MIN_f POWBLK_f PRODUCT PROD_f SIGNUM SINBLK_f SQRT SUMMATION SUM_f TANBLK_f TrigFun
" syn keyword scilabXcos			Matrix_pal CUMSUM EXTRACT EXTTRI MATBKSL MATCATH MATCATV MATDET MATDIAG MATDIV MATEIG MATEXPM MATINV
" syn keyword scilabXcos			MATLU MATMUL MATPINV MATRESH MATSING MATSUM MATTRAN MATZCONJ RICC ROOTCOEF SUBMAT
" syn keyword scilabXcos			Portaction_pal CLKINV_f CLKOUTV_f INIMPL_f IN_f OUTIMPL_f OUT_f
" syn keyword scilabXcos			Signalprocessing_pal QUANT_f SAMPHOLD_m
" syn keyword scilabXcos			Signalrouting_pal DEMUX EXTRACTOR FROM FROMMO GOTO GOTOMO GotoTagVisibility GotoTagVisibilityMO
" syn keyword scilabXcos			ISELECT_m MUX M_SWITCH NRMSOM_f RELAY_f SELECT_m SELF_SWITCH SWITCH2_m SWITCH_f
" syn keyword scilabXcos			Sinks_pal AFFICH_m BARXY CANIMXY CANIMXY3D CFSCOPE CMAT3D CMATVIEW CMSCOPE CSCOPE CSCOPXY CSCOPXY3D
" syn keyword scilabXcos			ENDBLK END_c TOWS_c TRASH_f WFILE_f WRITEAU_f WRITEC_f
" syn keyword scilabXcos			Sources_pal CLOCK_c CONST_m CURV_f Counter FROMWSB GENSIN_f GENSQR_f Modulo_Count PULSE_SC RAMP RAND_m
" syn keyword scilabXcos			READAU_f READC_f RFILE_f SAWTOOTH_f STEP_FUNCTION SampleCLK Sigbuilder TIME_f TKSCALE
" syn keyword scilabXcos			ThermoHydraulics_pal Bache Flowmeter PerteDP PuitsP SourceP VanneReglante
" syn keyword scilabXcos			Userdefinedfunctions_pal CBLOCK DSUPER EXPRESSION MBLOCK SUPER_f c_block fortran_block generic_block3 scifunc_block_m
" syn keyword scilabXcos			Zerocrossingdetection_pal GENERAL_f NEGTOPOS_f POSTONEG_f ZCROSS_f
" syn keyword scilabXcos			sci_struct
" syn keyword scilabXcos			curblock getblocklabel getscicosvars phase_simulation pointer_xproperty scicos_time set_blockerror set_xproperty
" syn keyword scilabXcos			scicos_block scicos_graphics scicos_model
" syn keyword scilabXcos			scicos_cpr scicos_sim scicos_state
" syn keyword scilabXcos			scicos_diagram scicos_params
" syn keyword scilabXcos			scicos_link
" syn keyword scilabXcos			block_parameter_error buildouttb create_palette getModelicaPath importXcosDiagram loadScicos
" syn keyword scilabXcos			loadXcosLibs scicos_debug scicos_getvalue standard_inputs standard_origin standard_outputs
" syn keyword scilabXcos			var2vec vec2var xcosPal xcosPalAdd xcosPalAddBlock xcosPalDelete xcosPalExport
" syn keyword scilabXcos			xcosPalGenerateAllIcons xcosPalMove
" syn keyword scilabXcos			xcos
syn keyword scilabSpreadsheet		csvDefault csvRead csvTextScan csvWrite read_csv readxls write_csv xls_open xls_read
syn keyword scilabConsole		clc completion diary lines prompt tohome
syn keyword scilabHistoryManager	addhistory browsehistory displayhistory gethistory gethistoryfile historymanager historysize loadhistory
syn keyword scilabHistoryManager	removelinehistory resethistory saveafterncommands saveconsecutivecommands savehistory sethistoryfile
syn keyword scilabMatlabBinaryFilesIO	loadmatfile matfile_close matfile_listvar matfile_open matfile_varreadnext matfile_varwrite savematfile
syn keyword scilabCompatibilityFunc	firstnonsingleton makecell mstr2sci mtlb_0 mtlb_a mtlb_all mtlb_any mtlb_axis mtlb_beta mtlb_box
syn keyword scilabCompatibilityFunc	mtlb_close mtlb_colordef mtlb_cumprod mtlb_cumsum mtlb_dec2hex mtlb_delete mtlb_diag mtlb_diff
syn keyword scilabCompatibilityFunc	mtlb_dir mtlb_double mtlb_e mtlb_echo mtlb_eval mtlb_exist mtlb_eye mtlb_false mtlb_fft
syn keyword scilabCompatibilityFunc	mtlb_fftshift mtlb_find mtlb_findstr mtlb_fliplr mtlb_fopen mtlb_format mtlb_fprintf mtlb_fread
syn keyword scilabCompatibilityFunc	mtlb_fscanf mtlb_full mtlb_fwrite mtlb_grid mtlb_hold mtlb_i mtlb_ifft mtlb_imp mtlb_int16
syn keyword scilabCompatibilityFunc	mtlb_int32 mtlb_int8 mtlb_is mtlb_isa mtlb_isfield mtlb_isletter mtlb_isspace mtlb_l
syn keyword scilabCompatibilityFunc	mtlb_legendre mtlb_linspace mtlb_logic mtlb_logical mtlb_lower mtlb_max mtlb_min mtlb_mode
syn keyword scilabCompatibilityFunc	mtlb_more mtlb_num2str mtlb_ones mtlb_plot mtlb_prod mtlb_rand mtlb_randn mtlb_rcond
syn keyword scilabCompatibilityFunc	mtlb_realmax mtlb_realmin mtlb_s mtlb_setstr mtlb_size mtlb_sort mtlb_strcmp mtlb_strcmpi
syn keyword scilabCompatibilityFunc	mtlb_strfind mtlb_strrep mtlb_sum mtlb_t mtlb_toeplitz mtlb_tril mtlb_triu mtlb_true
syn keyword scilabCompatibilityFunc	mtlb_uint16 mtlb_uint32 mtlb_uint8 mtlb_upper mtlb_var mtlb_zeros
syn keyword scilabAdvancedFunctions	clearfun funptr intppty newfun readgateway what
syn keyword scilabAdvancedFunctions	genlib get_function_path lib librarieslist libraryinfo whereis
syn keyword scilabAdvancedFunctions	add_profiling plotprofile profile remove_profiling reset_profiling showprofile
syn keyword scilabAdvancedFunctions	argn bytecode bytecodewalk code2str comp deff edit exec execstr fun2string funcprot
syn keyword scilabAdvancedFunctions	getd head_comments listfunctions macr2lst macr2tree macro macrovar mode recompilefunction
syn keyword scilabAdvancedFunctions	sciargs str2code tree2code
syn keyword scilabDevelopmentTools	assert_checkalmostequal assert_checkequal assert_checkerror assert_checkfalse assert_checkfilesequal
syn keyword scilabDevelopmentTools	assert_checktrue assert_comparecomplex assert_computedigits assert_cond2reltol assert_cond2reqdigits assert_generror
syn keyword scilabDevelopmentTools	bench_run example_run test_run user
syn keyword scilabDemoTools		add_demo demo_begin demo_choose demo_compiler demo_end demo_file_choice demo_function_choice demo_mdialog demo_message demo_run
syn keyword scilabDevelopmentTools	G_make addinter c_link call chooselcccompiler configure_ifort configure_msvc dllinfo
syn keyword scilabDevelopmentTools	findmsifortcompiler findmsvccompiler fort getdynlibext haveacompiler ilib_build
syn keyword scilabDevelopmentTools	ilib_compile ilib_for_link ilib_gen_Make ilib_gen_cleaner ilib_gen_gateway
syn keyword scilabDevelopmentTools	ilib_gen_loader ilib_include_flag ilib_mex_build ilib_verbose link ulink
syn keyword scilabATOMS			atomsAutoload atomsAutoloadAdd atomsAutoloadDel atomsAutoloadList atomsCategoryList atomsCheckModule
syn keyword scilabATOMS			atomsDepTreeShow atomsGetConfig atomsGetInstalled atomsGetInstalledPath atomsGetLoaded atomsGetLoadedPath
syn keyword scilabATOMS			atomsInstall atomsIsInstalled atomsIsLoaded atomsList atomsLoad atomsQuit atomsRemove atomsRepositoryAdd
syn keyword scilabATOMS			atomsRepositoryDel atomsRepositoryList atomsRestoreConfig atomsSaveConfig atomsSearch atomsSetConfig
syn keyword scilabATOMS			atomsShow atomsSystemInit atomsSystemUpdate atomsTest atomsUpdate atomsVersion
syn keyword scilabTclTKInterface	ScilabEval TCL_CreateSlave TCL_DeleteInterp TCL_EvalFile TCL_EvalStr TCL_ExistArray TCL_ExistInterp
syn keyword scilabTclTKInterface	TCL_ExistVar TCL_GetVar TCL_GetVersion TCL_SetVar TCL_UnsetVar TCL_UpVar winclose winlist
syn keyword scilabTextEditor		edit_error editor scinotes
syn keyword scilabUIData		browsevar closeEditvar editvar filebrowser
syn keyword scilabOnlineHelpManagement	add_help_chapter apropos del_help_chapter help help_from_sci help_skeleton manedit xmltochm xmltohtml xmltojar xmltopdf xmltops %helps
syn keyword scilabParallel		parallel_concurrency parallel_run
syn keyword scilabModulesManager	tbx_build_blocks tbx_build_cleaner tbx_build_gateway tbx_build_gateway_clean tbx_build_gateway_loader
syn keyword scilabModulesManager	tbx_build_help tbx_build_help_loader tbx_build_loader tbx_build_macros tbx_build_src tbx_builder_gateway
syn keyword scilabModulesManager	tbx_builder_gateway_lang tbx_builder_help tbx_builder_help_lang tbx_builder_macros tbx_builder_src tbx_builder_src_lang
syn keyword scilabLocalization		dgettext getdefaultlanguage getlanguage setdefaultlanguage setlanguage gettext
syn match scilabLocalization		"=\s*\zs_\ze(.*)"
syn keyword scilabJVM			javaclasspath javalibrarypath jre_path system_getproperty system_setproperty with_embedded_jre
syn keyword scilabAPI			isBooleanType getScalarBoolean createScalarBoolean
syn keyword scilabAPI			isBooleanSparseType getAllocatedBooleanSparseMatrix freeAllocatedBooleanSparse
syn keyword scilabAPI			CheckLhs (deprecated) CheckRhs (deprecated) Lhs (deprecated) LhsVar (deprecated) Rhs (deprecated) Scierror sci_types sciprint
syn keyword scilabAPI			isDoubleType getScalarDouble getScalarComplexDouble createScalarDouble createScalarComplexDouble
syn keyword scilabAPI			isIntegerType getScalarInteger8 createScalarInteger8
syn keyword scilabAPI			isListType isNamedListType isTListType isNamedTListType isMListType isNamedMListType
syn keyword scilabAPI			AssignOutputVariable CallOverloadFunction CheckInputArgument CheckOutputArgument ReturnArguments deleteNamedVariable
syn keyword scilabAPI			isPolyType getAllocatedSinglePoly getAllocatedSinglePoly getAllocatedMatrixOfPoly getAllocatedMatrixOfComplexPoly
syn keyword scilabAPI			freeAllocatedSinglePoly freeAllocatedSingleComplexPoly freeAllocatedSinglePoly freeAllocatedSinglePoly
syn keyword scilabAPI			isSparseType getAllocatedSparseMatrix getAllocatedComplexSparseMatrix freeAllocatedSparseMatrix freeAllocatedComplexSparseMatrix
syn keyword scilabAPI			isStringType getAllocatedSingleString getAllocatedMatrixOfString createSingleString freeAllocatedSingleString freeAllocatedMatrixOfString
syn keyword scilabAPI			api_scilab
syn keyword scilabCallScilabAPI		DisableInteractiveMode GetLastJob ScilabHaveAGraph SendScilabJob SendScilabJobs StartScilab TerminateScilab call_scilab fromc fromjava
syn keyword scilabPreferences		preferences
syn keyword scilabWindowsTools		consolebox createGUID dos findfileassociation getsystemmetrics istssession mcisendstring powershell win64 winopen winqueryreg

" Define the default highlighting.
" For version 5.7 and earlier: only when not done already
" For version 5.8 and later: only when an item doesn't have highlighting yet
if version >= 508 || !exists("did_scilab_syntax_inits")
  if version < 508
    let did_scilab_syntax_inits = 1
    command -nargs=+ HiLink hi link <args>
  else
    command -nargs=+ HiLink hi def link <args>
  endif

  HiLink scilabOperator			Operator
  HiLink scilabTransposeOperator	Operator
  HiLink scilabSizeOperator		Float
  HiLink scilabFloat			Float
  HiLink scilabAPI			Function
  HiLink scilabARnoldiPACKage		Function
  HiLink scilabATOMS			Function
  HiLink scilabAdvancedFunctions	Function
  HiLink scilabBoolean			Function
  HiLink scilabCACSD			Function
  HiLink scilabCallScilabAPI		Function
  HiLink scilabCompatibilityFunc	Function
  HiLink scilabConsole			Function
  HiLink scilabDataStructure		Function
  HiLink scilabDemoTools		Function
  HiLink scilabDevelopmentTools		Function
  HiLink scilabDiffAndIntegral		Function
  HiLink scilabElementaryFunction	Function
  HiLink scilabFFTW			Function
  HiLink scilabFilesIO			Function
  HiLink scilabGUI			Function
  HiLink scilabGeneticAlgorithms	Function
  HiLink scilabGraphics			Function
  HiLink scilabHistoryManager		Function
  HiLink scilabIOFunctions		Function
  HiLink scilabIntegers			Function
  HiLink scilabInterpolation		Function
  HiLink scilabJVM			Function
  HiLink scilabLinearAlgebra		Function
  HiLink scilabLocalization		Function
  HiLink scilabMatlabBinaryFilesIO	Function
  HiLink scilabModulesManager		Function
  HiLink scilabOnlineHelpManagement	Function
  HiLink scilabOptAndSim		Function
  HiLink scilabOutputFunctions		Function
  HiLink scilabParallel			Function
  HiLink scilabParameters		Function
  HiLink scilabPolynomials		Function
  HiLink scilabPreferences		Function
  HiLink scilabRandlib			Function
  HiLink scilabScilab			Function
  HiLink scilabSignalProcessing		Function
  HiLink scilabSimulatedAnnealing	Function
  HiLink scilabSoundFileHandling	Function
  HiLink scilabSparseMatrix		Function
  HiLink scilabSpecialFunctions		Function
  HiLink scilabSpreadsheet		Function
  HiLink scilabStatistics		Function
  HiLink scilabStrings			Function
  HiLink scilabTclTKInterface		Function
  HiLink scilabTextEditor		Function
  HiLink scilabTimeAndDate		Function
  HiLink scilabUIData			Function
  HiLink scilabUMFPACKInterface		Function
  HiLink scilabWindowsTools		Function
  HiLink scilabXMLManagement		Function
"   HiLink scilabXcos			Function
  HiLink scilabLineContinuation		Special
  HiLink scilabLabel			Label
  HiLink scilabConditional		Conditional
  HiLink scilabRepeat			Repeat
  HiLink scilabString			String
  HiLink scilabDelimiter		Delimiter
  HiLink scilabTransposeOther		Identifier
  HiLink scilabSpecialVariables		Identifier
  HiLink scilabNumber			Number
  HiLink scilabError			Error
  HiLink scilabTodo			Todo
  HiLink scilabStatement		Statement
  HiLink scilabSemicolon		SpecialChar
  HiLink scilabColon			SpecialChar
  HiLink scilabComment			Comment
  HiLink scilabArithmeticOperator	scilabOperator
  HiLink scilabRelationalOperator	scilabOperator
  HiLink scilabLogicalOperator		scilabOperator
  HiLink scilabPreDefVariables		Constant

"optional highlighting
  "HiLink scilabIdentifier		Identifier
  "HiLink scilabTab			Error

  delcommand HiLink
endif

let b:current_syntax = "scilab"

"EOF	vim: ts=8 noet tw=100 sw=8 sts=0
