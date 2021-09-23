#include "CustomMath.h"

#include <cmath>
#include <cstdint>
#include <cstring>


const double
customMath::taylorTableLn[64][8] =  {{0.70092932100200022738, 0.4961240310077519311, -0.12306952707169040162, 0.040705166576683132762, -0.01514610849364953632, 0.006011478719960126893, -0.002485365879053282219, 0.0010568997758100336527},
                                      {0.71631423984147968298, 0.48854961832061066795, -0.11934036478060719433, 0.038869126442538724786, -0.014242122665968389231, 0.0055663868740273396429, -0.0022662134856090950311, 0.00094899234293772353078},
                                      {0.73146604486208188778, 0.48120300751879696577, -0.11577816722256770166, 0.03714186818167584947, -0.013404584005416849898, 0.0051602609103559606007, -0.0020692775580374777808, 0.00085349364370181665749},
                                      {0.7463916950787575777, 0.47407407407407409217, -0.11237311385459533564, 0.035515453267625188283, -0.012627716717377846997, 0.0047891784883684869814, -0.0018920211312073035095, 0.00076882128506201537817},
                                      {0.76109784246845302302, 0.46715328467153283132, -0.10911609568970109807, 0.033982628341318099618, -0.011906322338564007124, 0.0044496620710545780816, -0.001732228543476234655, 0.00069361393190289267769},
                                      {0.77559084977101988567, 0.46043165467625901677, -0.10599865431395889825, 0.032536757199488577919, -0.01123571471637015845, 0.0041386229746629655743, -0.0015879608535877081342, 0.00062669780861015404825},
                                      {0.7898768070184963852, 0.453900709219858145, -0.10301292691514511712, 0.031171760390398518814, -0.01061166311162502901,0.0038533131298950462455, -0.0014575179687546272692, 0.00056705866261578207182},
                                      {0.80396154690023546863, 0.44755244755244755206, -0.10015159665509315579, 0.029882061472848307077, -0.010030342312564466747, 0.0035912834014216834556, -0.0013394063968006277455, 0.00051381823813330775638},
                                      {0.81785065906090259613, 0.44137931034482758008, -0.097407847800237812308, 0.028662539122828595284, -0.0094882888130742946281, 0.0033503474981338198305, -0.0012323117234515198679, 0.00046621448453732373517},
                                      {0.83154950341906441746, 0.43537414965986392934, -0.094775325096024806704, 0.027508484381612637037, -0.008982362247057189622, 0.0031285506601995112866, -0.0011350750694601398263, 0.00042358486557113089917},
                                      {0.84506322258578725481, 0.42953020134228186988, -0.092248096932570605722, 0.026415562432592923359, -0.0085097113876809425681, 0.0029241424365722434343, -0.0010466729079453220182, 0.00038535224990508495403},
                                      {0.8583967534552524592, 0.42384105960264900625, -0.089820621902548136717, 0.025379778374229933907, -0.0080677441189605098026, 0.0027355529727866099814, -0.00096619972548754433101, 0.00035101295609008231454},
                                      {0.87155483803276356802, 0.41830065359477125453, -0.087487718398906400008, 0.024397446525185224886, -0.0076541008706463463038, 0.0025613723174973397449, -0.00089285309542390910165, 0.00032012660003994498166},
                                      {0.88454203355957472521, 0.41290322580645161255, -0.085244536940686788107, 0.023465162856791198231, -0.007266631078232113454, 0.0024003323303579630205, -0.0008259208018436001511, 0.00029230745429303451113},
                                      {0.89736272198863620275, 0.40764331210191084853, -0.083086534950707932601, 0.022579780198918503881, -0.0069033722901152118692, 0.0022512908360120945962, -0.00076476971074720837083, 0.00026721707818646770307},
                                      {0.91002111886055969681, 0.40251572327044027322, -0.081009453739962822283, 0.021738385909256269929, -0.0065625315952471761222, 0.002113217721236826472, -0.00070883613290962732903, 0.00024455801890143475091},
                                      {0.92252128162479118956, 0.39751552795031053211, -0.079009297480807066938, 0.020938281734044109778, -0.0062424690884106664307, 0.0019851827163144481272, -0.0006576174629613078422, 0.00022406841683863547754},
                                      {0.93486711744709050098, 0.39263803680981596012, -0.07708231397493318926, 0.020176965621250403382, -0.0059416831277301811534, 0.0018663446388943881644, -0.00061066491252168117639, 0.00020551737634384362576},
                                      {0.94706239054090879392, 0.38787878787878787845, -0.075224977043158863799, 0.019452115275806734079, -0.0056587971711437781014, 0.00175594191007612992, -0.00056757718305491062973, 0.00018870098553513909917},
                                      {0.95911072905708327774, 0.3832335329341317598, -0.073433970382588115156, 0.018761573271399756974, -0.0053925480061508285906, 0.001653284179131271972, -0.00052799494742715273599, 0.00017343888777761047419},
                                      {0.97101563156340164884, 0.37869822485207099705, -0.071706172753054867042, 0.018103333555011875966, -0.0051417752108909482722, 0.0015577449159622281695, -0.00049159602870602850822, 0.00015957132292740064785},
                                      {0.98278047314298799808, 0.3742690058479531956, -0.070038644369207619933, 0.017475529199334455255, -0.0049054117050763386063, 0.0014687548497070673607, -0.00045809118119518666464, 0.00014695656940597466209},
                                      {0.9944085111381071318, 0.36994219653179188922, -0.06842861438738347124, 0.01687642127473041459, -0.004682475266977224361, 0.0013857961483770746356, -0.00042722039256325996734, 0.00013546872893830869187},
                                      {1.0059028905638423002, 0.36571428571428571397, -0.066873469387755096749, 0.016304388726919340213, -0.0044720609079550183856, 0.0013083972484988397669, -0.00039874963763774159834, 0.00012499580477787165453},
                                      {1.0172666492141573258, 0.36158192090395480101, -0.065370742762296915607, 0.015757919159273076071, -0.0042733340092943944075, 0.0012361282557958926803, -0.00037246802434527838269, 0.00011543803175834291893},
                                      {1.0285027224810832092, 0.3575418994413407936, -0.063918104928060923142, 0.015235600429779883741, -0.0040855241375946059479, 0.0011685968482952168845, -0.000348185280683863509, 0.00010670642281133564384},
                                      {1.0396139479061539124, 0.35359116022099446042, -0.062513354293214487556, 0.014736112982562530466, -0.003907919464988958945, 0.001105444622140523258, -0.00032572953875227204056, 9.8721501879141631262e-05},
                                      {1.0506030694817489746, 0.3497267759562841527, -0.061154408910388487197, 0.014258222842494946286, -0.0037398617291790024753, 0.0010463438280544532736, -0.00030494537793936701917, 9.1412197602433207318e-05},
                                      {1.0614727417186529124, 0.34594594594594596515, -0.059839298758217675245, 0.013800775209102453603, -0.0035807416758752317817, 0.00099099445299898297583, -0.00028569209455826533994, 8.4714875915346644973e-05},
                                      {1.0722255334949146732, 0.34224598930481281434, -0.058566158597615027015, 0.013362688592682216843, -0.0034299949328809968552, 0.00093912160729148165192, -0.00026784216963749915823, 7.8572492850114349864e-05},
                                      {1.0828639316999704167, 0.33862433862433860554, -0.057333221354385378865, 0.01294294944155436819, -0.0032870982708709511194, 0.00089047318237350640314, -0.00025127991036994889052, 7.2933851535948878331e-05},
                                      {1.093390344686957949, 0.33507853403141363291, -0.056138811984320607928, 0.012540607214647534612, -0.0031515662110108986008, 0.00084481774871077500427, -0.00023590024396810468509, 6.7752949651273143951e-05},
                                      {1.1038071055452136626, 0.33160621761658032325, -0.054981341780987409706, 0.01215476985831846074, -0.0030229479440377518377, 0.00080194266701934149421, -0.00022160764546303043064, 6.2988405520209982428e-05},
                                      {1.1141164752040748631, 0.32820512820512820484, -0.053859303090072321862, 0.01178459965047736098, -0.0029008245293482737455, 0.00076165238924426478199, -0.00020831518338304674462, 5.860295268797798457e-05},
                                      {1.1243206453783165522, 0.32487309644670048225, -0.052771264397433584326, 0.011429309378801183472, -0.0027848063461038412465, 0.00072376692853054155376, -0.00019594366931452900267, 5.4562994210862322743e-05},
                                      {1.1344217413648205461, 0.32160804020100502987, -0.051715865760965630538, 0.011088158823121608068, -0.002674530771406217599, 0.00068812047987938868617, -0.0001844208991301543902, 5.083820909259100971e-05},
                                      {1.1444218246994037891, 0.31840796019900496905, -0.050691814559045568489, 0.010760451515021281568, -0.0025696600632886643263, 0.00065456017532527176455, -0.00017368097521068570112, 4.7401204321892896566e-05},
                                      {1.1543228956821154352, 0.31527093596059113656, -0.049697881530733574451, 0.010445531750301967344, -0.0024698794286428294754, 0.00062294495934242803854, -0.00016366370031984972706, 4.422720684223946456e-05},
                                      {1.1641268957787362925, 0.31219512195121951192, -0.048732897085068414833, 0.010142781832339440445, -0.0023748952583038693102, 0.00059314457183003957592, -0.00015431403494765256004, 4.1293790536514694139e-05},
                                      {1.1738357099056972999, 0.30917874396135264226, -0.047795747858759833615, 0.0098516195264432490641, -0.002284433513378148928, 0.00056503862746357981246, -0.00014558161094069689123, 3.8580633955298555403e-05},
                                      {1.1834511686051392143, 0.30622009569377989235, -0.046885373503353859614, 0.0095714957072237545621, -0.0021982382485489964166, 0.00053851578146272064112, -0.0001374202951101041281, 3.6069305073328771242e-05},
                                      {1.1929750501163947085, 0.30331753554502371983, -0.046000763684553358512, 0.0093018921826584985424, -0.0021160702595621229405, 0.00051347297293640142661, -0.00012978779726670493493, 3.3743069837789230019e-05},
                                      {1.2024090823497532998, 0.30046948356807512415, -0.045140955277832878201, 0.0090423196800666801043, -0.0020377058433953083513, 0.00048981473794290977808, -0.00012264531779478178466, 3.1586721685577602884e-05},
                                      {1.2117549447679909491, 0.29767441860465115866, -0.044305029745808542641, 0.0087923159805635558994, -0.0019629356607769798612, 0.00046745258526409942798, -0.0001159572304531099248, 2.9586429564115754949e-05},
                                      {1.2210142701807877863, 0.2949308755760368661, -0.043492110684023872758, 0.0085514441897927111202, -0.0018915636917513831399, 0.00044630442865286094282, -0.00010969079659670929279, 2.7729602299628945249e-05},
                                      {1.2301886464568290158, 0.292237442922374413, -0.042701361522904025814, 0.0083192911338382264513, -0.0018234062759097484623, 0.00042629406998437958725, -0.0001038159074543390021, 2.6004767424961628905e-05},
                                      {1.2392796181580809645, 0.28959276018099550098, -0.041931983374623782446, 0.0080954658702139421778, -0.0017582912297297249773, 0.00040735072833557431403, -9.8304851483698178684e-05, 2.4401462811725984264e-05},
                                      {1.2482886881004469082, 0.28699551569506726034, -0.041183213014538797014, 0.0078795983047249117809, -0.0016960570342008779703, 0.00038940861054298189, -9.3132104165735568884e-05, 2.291013965383885653e-05},
                                      {1.2572173188447481884, 0.2844444444444444442, -0.040454320987654321384, 0.0076713379058070416475, -0.0016365520865721688067, 0.00037240651925553357529, -8.8274137897607953267e-05, 2.1522075525512034165e-05},
                                      {1.2660669341217307693, 0.28193832599118945348, -0.039744609831357101404, 0.0074703525086838593397, -0.0015796340106468073367, 0.00035628749491240771291, -8.3709249905998856672e-05, 2.022929639012181207e-05},
                                      {1.2748389201945677929, 0.27947598253275107716, -0.039053412406323297079, 0.0072763272023426226603, -0.0015251690205783663152, 0.00034099848844372211785, -7.9417406333326823892e-05, 1.9024506570179354606e-05},
                                      {1.283534627162121744, 0.27705627705627705604, -0.038380090328142275891, 0.0070889632929324840338, -0.0014730313335963602511, 0.00032649006181876042386, -7.5380100852816251224e-05, 1.7901025805492539193e-05},
                                      {1.2921553702060286639, 0.27467811158798283167, -0.03772403249277017645, 0.0069079773377318772279, -0.001423102627515579954, 0.00031271611385750083941, -7.1580226347925797166e-05, 1.6852732628818824852e-05},
                                      {1.3007024307844869959, 0.27234042553191489811, -0.037084653689452240499, 0.0067331002443260798765, -0.0013752715392666037698, 0.0002996336289806387815, -6.8001958350215877718e-05, 1.5874013377801156485e-05},
                                      {1.3091770577754593408, 0.2700421940928269815, -0.03646139329523402639, 0.0065640764300843252102, -0.0013294332010297367844, 0.00028720244680473640351, -6.4630649069982869953e-05, 1.495971624043002988e-05},
                                      {1.3175804685718388587, 0.26778242677824265483, -0.035853714045622452067, 0.0064006630374332969258, -0.0012854888108652646854, 0.00027538505069582238424, -6.145273097814724127e-05, 1.4105109800124651735e-05},
                                      {1.3259138501309830627, 0.26556016597510373467, -0.035261100876362319601, 0.0062426292007944361728, -0.0012433452350129995962, 0.00026414637357952526703, -5.8455629007917063651e-05, 1.3305845607018465843e-05},
                                      {1.334178359980876527, 0.26337448559670784132, -0.034683059831665224992, 0.0060897553613897788607, -0.0012029146392868699921, 0.00025345361947114300244, -5.5627680542363342538e-05, 1.2557924355242519102e-05},
                                      {1.3423751271850550282, 0.26122448979591839091, -0.03411911703456892847, 0.005941832626428330405, -0.0011641141472186118146, 0.0002432760993371139806, -5.2958062440732289289e-05, 1.1857665292852009984e-05},
                                      {1.3505052532683052835, 0.25910931174089069096, -0.033568817715419037306, 0.0057986621694650960362, -0.0011268655228110309105, 0.00023358508003208414377, -5.0436724433378495253e-05, 1.1201678532340855684e-05},
                                      {1.3585698131050356618, 0.25702811244979917316, -0.033031725294753309929, 0.0056600546689804867584, -0.0010910948759480455696, 0.00022435364517485922777, -4.8054328283771712449e-05, 1.0586839966132149355e-05},
                                      {1.3665698557721119855, 0.25498007968127489598, -0.032507420517134649751, 0.0055258297824611352783, -0.0010567323886778267998, 0.000215556566933484961, -4.580219217710171923e-05, 1.0010268523623825325e-05},
                                      {1.3745064053678484051, 0.25296442687747033862, -0.031995500632723522516, 0.0053958156534764309398, -0.0010237120607386113787, 0.00020717018778583756299, -4.3672239849451915884e-05, 9.4693055348331654037e-06},
                                      {1.3823804617987542898, 0.25098039215686274161, -0.031495578623606308721, 0.0052698484494400083672, -0.00099197147283576654619, 0.00019917231140859312255, -4.1656954020097901694e-05, 8.9614959908782034897e-06}};

const double
customMath::taylorWeightsCos[] = {0.16666666666666665741, -0.0083333333333333332177, 0.00019841269841269841253,
                                  -2.7557319223985892511e-06, 2.5052108385441720224e-08, -1.6059043836821613341e-10,
                                  7.6471637318198164055e-13, -2.8114572543455205981e-15, 8.2206352466243294955e-18};


double customMath::ln(double x)
{
    uint64_t xBits;
    double x1;
    int exponent;
    // basic principle, ln(a) = ln(b * 2^n) = ln(b) + n * ln(2), b element [2.0, 4.0], a is positive
    // converting double to an int so bitwise operations can be performed
    std::memcpy(&xBits, &x, sizeof(x));
    // extracting the exponent bits and manipulating them to their final form
    unsigned long long exp = xBits & 0x7ff0000000000000;
    exp >>= 52;
    exp -= 1024;
    exponent = (int) exp;
    // changing the exponent so all numbers land in [2.0, 4.0]
    xBits &= 0x400fffffffffffff;
    xBits |= 0x4000000000000000;
    // converting the altered number back to double
    std::memcpy(&x1, &xBits, sizeof(xBits));

    double output = exponent * 0.69314718055994528623;

    // faster version
    int row = (int) ((x1 - 2.0) * 32.0);
    auto usedRow = taylorTableLn[row];
    x1 -= (2 + (row + 0.5) * 0.03125);

    double x2 = x1 * x1;
    double x4 = x2 * x2;
    double x6 = x4 * x2;

    output += usedRow[0] + usedRow[1] * x1;
    output += usedRow[2] * x2 + usedRow[3] * x2 * x1;
    output += usedRow[4] * x4 + usedRow[5] * x4 * x1;
    output += usedRow[6] * x6 + usedRow[7] * x6 * x1;

    return output;
}

double customMath::cos(double x)
{
    double x1 = x - M_PI_2;
    double output = 0.0;

    double x2 = x1 * x1;
    double x3 = x1 * x2;
    double x5 = x3 * x2;
    double x7 = x5 * x2;
    double x9 = x7 * x2;
    double x11 = x9 * x2;
    double x13 = x11 * x2;
    double x15 = x13 * x2;
    double x17 = x15 * x2;
    double x19 = x17 * x2;

    output += -x1 + taylorWeightsCos[0] * x3;
    output += taylorWeightsCos[1] * x5 + taylorWeightsCos[2] * x7;
    output += taylorWeightsCos[3] * x9 + taylorWeightsCos[4] * x11;
    output += taylorWeightsCos[5] * x13 + taylorWeightsCos[6] * x15;
    output += taylorWeightsCos[7] * x17 + taylorWeightsCos[8] * x19;

    return output;
}

float customMath::lnf(float x)
{
    auto xd = (double) x;
    uint64_t xBits;
    double x1;
    int exponent;
    // basic principle, ln(a) = ln(b * 2^n) = ln(b) + n * ln(2), b element [2.0, 4.0], a is positive
    // converting double to an int so bitwise operations can be performed
    std::memcpy(&xBits, &xd, sizeof(xd));
    // extracting the exponent bits and manipulating them to their final form
    unsigned long long exp = xBits & 0x7ff0000000000000;
    exp >>= 52;
    exp -= 1024;
    exponent = (int) exp;
    // changing the exponent so all numbers land in [2.0, 4.0]
    xBits &= 0x400fffffffffffff;
    xBits |= 0x4000000000000000;
    // converting the altered number back to double
    std::memcpy(&x1, &xBits, sizeof(xBits));

    double output = exponent * 0.69314718055994528623;
    int row = (int) ((x1 - 2.0) * 32.0);
    x1 -= (2 + (row + 0.5) * 0.03125);

    double x2 = x1 * x1;

    output += taylorTableLn[row][0] + taylorTableLn[row][1] * x1;
    output += taylorTableLn[row][2] * x2 + taylorTableLn[row][3] * x2 * x1;

    return (float) output;
}

float customMath::cosf(float x)
{
    double x1 = x - M_PI_2;
    double output = 0.0;

    double x2 = x1 * x1;
    double x3 = x1 * x2;
    double x5 = x3 * x2;
    double x7 = x5 * x2;
    double x9 = x7 * x2;

    output += -x1 + taylorWeightsCos[0] * x3;
    output += taylorWeightsCos[1] * x5 + taylorWeightsCos[2] * x7;
    output += taylorWeightsCos[3] * x9 + taylorWeightsCos[4] * x9 * x2;

    return (float) output;
}