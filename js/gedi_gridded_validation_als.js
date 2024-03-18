/* 
Title: gedi_gridded_validation_als

By: Patrick Burns [pb463@nau.edu], Northern Arizona University

About: validate gridded GEDI rasters with ALS 

*/



// ----------------------
// ----- IMPORTS --------
// ----------------------
var palettes = require('users/gena/packages:palettes');

var borneo = 
    /* color: #2d8ad6 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[109.31324050561682, -0.9513194344818261],
          [109.36542556421057, -1.3522385721526273],
          [109.90235757420857, -1.3110454997308683],
          [110.0328437906502, -1.4318498984930932],
          [109.76917312280432, -1.8546646594925995],
          [109.95594070092932, -1.9946612847236098],
          [110.08296028594235, -2.376783457033743],
          [110.15438125268716, -2.6992502365184463],
          [110.21102502595278, -2.885434253298825],
          [110.23609206811685, -2.994143127733226],
          [110.48877761499185, -3.0277423100492156],
          [111.16169021264807, -3.138817839755338],
          [111.59427688256997, -3.0489983981148088],
          [111.66706130639807, -3.2923842198590845],
          [111.77692458764807, -3.6021860195347477],
          [112.17243240014807, -3.6076683144090205],
          [113.06369826928872, -3.448668538297063],
          [114.06070754663247, -3.5048698366824413],
          [114.44110915796057, -3.659748459117849],
          [114.53449294702307, -4.229681705281973],
          [115.06458327905435, -4.137916411676063],
          [115.90778396264807, -3.7364926150497157],
          [116.01007902402809, -3.8406126013966886],
          [116.03275344506997, -4.139286121680654],
          [116.19551822981528, -4.016719734193729],
          [116.38019607202307, -4.070797734787796],
          [116.5065388454606, -3.148416411278458],
          [116.44338131891085, -2.8034604909142087],
          [116.61090896264807, -2.5353132071949815],
          [116.61090896264807, -1.8985864012334772],
          [116.53950348148435, -1.7338669362613965],
          [116.73175857202307, -1.5416922149140027],
          [117.69306228296057, -0.8387233284552714],
          [117.73700759546057, -0.38281781891186184],
          [117.64362380639807, -0.047737653869886834],
          [117.82215163842935, 0.5922054092376124],
          [118.00204501322762, 0.756314855801564],
          [118.25610164156642, 0.7199296612565214],
          [118.79718825952307, 0.6938225314359999],
          [119.09931228296057, 0.9574667979091466],
          [118.8685993923356, 1.2265827022085327],
          [118.60429052447193, 1.4991083799084075],
          [118.0940632595231, 1.7633547541950694],
          [117.96635619508878, 1.9486249332455403],
          [118.12429789766853, 2.2244785392276007],
          [118.10368811186878, 2.369247870883394],
          [118.11054275171057, 2.5057844538253],
          [117.89084895357006, 2.6402178236861658],
          [117.78644607202307, 2.988625898545005],
          [117.73148998843229, 3.1709691801121225],
          [117.71362623957441, 3.591154107459863],
          [117.84549758569497, 3.6549202006824824],
          [117.92377517358557, 3.9043153449856827],
          [117.98763320581216, 4.135145473606968],
          [118.28701064721841, 4.228965422436078],
          [118.53419238342074, 4.332361981598213],
          [118.63512527494152, 4.404260650455562],
          [118.7752156032731, 4.443270168089994],
          [118.58374395108005, 4.674316970645196],
          [118.47144465043885, 4.73128507674661],
          [118.38191141426582, 4.729104946073068],
          [118.28550825730417, 4.75702160255018],
          [118.24924514428872, 4.82653348181836],
          [119.24213454858557, 4.988671818222113],
          [119.41791579858557, 5.2622315307124765],
          [119.30255935327307, 5.606744804159113],
          [118.94005765125252, 5.500196486489413],
          [118.66397903100747, 5.698307893240144],
          [118.37970779077307, 5.951054837597885],
          [118.00484824333341, 6.0958480824518695],
          [117.93615867668454, 6.036103653428742],
          [117.74730727807778, 5.985883929768664],
          [117.76172683374185, 6.469841000674811],
          [117.57018853539785, 6.592609002270759],
          [117.47264907495278, 6.774723330162043],
          [117.38339193799064, 6.774724384605027],
          [117.2728405051289, 6.771312013847622],
          [117.25172331592925, 6.7881870737684915],
          [117.24502441073479, 6.815289163605551],
          [117.2438215082243, 6.8400011337778555],
          [117.26107319890983, 6.8738324446271095],
          [117.27764175073403, 6.938339614563771],
          [117.19765724626748, 6.968695109486998],
          [117.17997644623054, 6.999204347551059],
          [117.1519931123315, 7.0269863673354775],
          [116.88625381128091, 7.033074804620795],
          [116.67957351342935, 7.050792918701677],
          [116.5118776598698, 6.649803924606485],
          [116.14948318139807, 6.333368207010038],
          [116.08434868504074, 6.121716870589651],
          [116.08269166475871, 6.058899996929677],
          [116.07554248867086, 5.998131515346811],
          [116.03288827413496, 5.955802507672481],
          [116.01633539946661, 5.833565150948698],
          [115.70797011987466, 5.711972779527299],
          [115.34130141870278, 5.496029967212734],
          [115.335682073675, 5.228728956077563],
          [115.34381327397912, 5.008509095001121],
          [115.17058965041176, 5.070769943661493],
          [114.84211013452307, 5.098109736814388],
          [114.46308181421057, 4.840901753482645],
          [114.30720497106665, 4.643843350052354],
          [113.96456886266824, 4.602094242877821],
          [113.93778968786353, 4.490523198797702],
          [113.8704984280979, 4.319368952363201],
          [113.5628812405979, 3.987905904027598],
          [113.3431546780979, 3.754979334889153],
          [113.23603797887915, 3.5713343897969008],
          [112.92842079137915, 3.30539575499658],
          [112.35713172887915, 3.1353766440792374],
          [111.6732328030979, 2.959844070355028],
          [111.2008206937229, 2.817203301600053],
          [111.19807411169165, 2.6635706461874746],
          [111.03327918981665, 2.1092444115266105],
          [110.95291138144285, 1.5732501352993764],
          [110.59382606481668, 1.6247314693276416],
          [110.50799790235672, 1.7816076792304507],
          [110.33701633148837, 1.816266517128556],
          [110.0788419339573, 1.7194481903436285],
          [109.76987738026853, 1.8937672652595903],
          [109.67646766637915, 2.1559039659281987],
          [109.33039833044168, 1.9829822179331984],
          [109.21435523962137, 1.7455287919390905],
          [108.9925687405979, 1.5080453724376481],
          [108.8881986234104, 1.2554335048602767],
          [108.8112943265354, 0.9890663918247684],
          [108.81953407262915, 0.6622575093635708],
          [108.87240612719256, 0.4878562434357529],
          [108.84699989294165, 0.3216943941210404],
          [109.06947303747293, 0.13630166864840107],
          [109.0529935452854, -0.013386923372405095],
          [108.9925687405979, -0.1919143966633088],
          [108.9870755765354, -0.39790520935982016],
          [109.0694730374729, -0.6341015178215467],
          [109.1628568265354, -0.6807903372186308],
          [109.1628568265354, -0.8730333144914835]]]),
    test = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-111.73512386909013, 35.126834140277474],
          [-111.73512386909013, 35.11503981552214],
          [-111.72139095893388, 35.11503981552214],
          [-111.72139095893388, 35.126834140277474]]], null, false),
    jambi_als_geo = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[102.55796016698241, -1.6480784198342942],
          [102.55796016698241, -2.2656959759823114],
          [103.41214717870116, -2.2656959759823114],
          [103.41214717870116, -1.6480784198342942]]], null, false),
    borneo_safe_geo = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[116.90082837537796, 5.097031371801567],
          [116.90082837537796, 4.445613680241448],
          [117.9472761292842, 4.445613680241448],
          [117.9472761292842, 5.097031371801567]]], null, false);


// ----------------------
// ----- INPUTS ---------
// ----------------------
// ancillary datasets
// elevation for slope
var usned_dem = ee.Image("USGS/3DEP/10m")
var glo30_proj = ee.ImageCollection("COPERNICUS/DEM/GLO30").first().select(0).projection();
var glo30_dem = ee.ImageCollection("COPERNICUS/DEM/GLO30").select('DEM').mosaic().setDefaultProjection(glo30_proj)

// masks
var nlcd_not_water_mask = ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD").select('landcover').mosaic().neq(11)
var nlcd_not_urban_mask = ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD").select('landcover').mosaic().lt(21).or(
  ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD").select('landcover').mosaic().gt(24))
var glad_not_water_mask = ee.ImageCollection('projects/glad/water/annual').filter(ee.Filter.gte('year', 2019))
                                                                       .filter(ee.Filter.lte('year', 2023))
                                                                       .mean()
                                                                       .lt(10)
var cop_not_urban_mask = ee.Image(ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global")
                           .filterDate('2019-01-01', '2019-12-31')
                           .first())
                           .select('discrete_classification')
                           .neq(50)
Map.addLayer(cop_not_urban_mask)         

// ALS CHMs (height should be in meters, convert here if necessary)
// 10m CHMs using maxpool with 12m radius
// source: https://zenodo.org/records/7885699
var spain_chm = ee.Image('users/pb463/als_spain_chm_10m') // acq started 2015?
var switz_chm = ee.Image('users/pb463/als_switz_chm_10m') // acq 2017-2021
var gabon_chm = ee.Image('users/pb463/als_gabon_lvis_chm_10m') // acq 2016
// var suma_chm = ee.Image('users/pb463/als_sumatra2020_chm_10m') // acq 2020


// finer res. (<= 3m) CHMs (height should be in meters, convert here if necessary)
// Borneo SAFE 
var borneo_mal_safe_chm = ee.Image('users/pb463/Borneo_Maliau_ALS_SAFE_CH_1m')
var borneo_dan_safe_chm = ee.Image('users/pb463/Borneo_Danum_ALS_SAFE_CH_1m') // acq 2014
var borneo_safe_chm = ee.Image('users/pb463/Borneo_SAFE_ALS_SAFE_CH_1m') // acq 2014
var borneo_safe_mos_chm = ee.ImageCollection([borneo_mal_safe_chm,
                                          borneo_dan_safe_chm,
                                          borneo_safe_chm])
                        .mosaic()
                        .clip(borneo_safe_geo)
                        .setDefaultProjection(borneo_safe_chm.projection()) // acq 2014

// other ALS metrics
// SAFE, https://zenodo.org/records/4020697
var borneo_pai_25m = ee.Image('users/pb463/als_safemos_pai_25m')
var borneo_padshan_25m = ee.Image('users/pb463/als_safemos_padshan_25m')

// NASA CMS Borneo   
var borneo_cms_chm = ee.Image('users/pb463/als_borneo_2014_chm_mos') // acq 2014

// Sonoma County
var sonoma_chm = ee.Image('users/pb463/als_ca-sonoma2013_chm_mos_3m') // acq 2013

// Coconino NF
var coco_chm = ee.Image("users/pb463/als_az-coco2019_chm_1m").divide(100) // acq 2019

// DFG Sumatra
var suma2020_chm = ee.Image('users/pb463/als_sumatradfg_2020_chm_mos')
var suma2022_chm = ee.Image('users/pb463/als_sumatradfg_2022_chm_mos')
var suma_chm = ee.ImageCollection([suma2020_chm, suma2022_chm])
                 .reduce(ee.Reducer.firstNonNull())
                 .setDefaultProjection(suma2020_chm.projection())
                 .clip(jambi_als_geo) // use acq 2020 as main year



// Jambi
var suma_rh50_25m = ee.Image("projects/earthengine-legacy/assets/users/pb463/als_suma_zq50_25m")
var suma_pai_25m = ee.Image("projects/earthengine-legacy/assets/users/pb463/als_suma_lai_25m")
var suma_fhd_25m = ee.Image("projects/earthengine-legacy/assets/users/pb463/als_suma_fhd_25m")



// GEDI gridded
//var gedi_rh95_1k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-95-a0_vf_20190417_20230316_1000m")
//var gedi_rh95_6k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-95-a0_vf_20190417_20230316_6000m")
//var gedi_rh95_12k = ee.Image('projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-95-a0_vf_20190417_20230316_12000m')
var gedi_rh50_1k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-50-a0_vf_20190417_20230316_1000m")
// var gedi_rh50_6k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-50-a0_vf_20190417_20230316_6000m")
var gedi_rh98_1k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-98-a0_vf_20190417_20230316_1000m")
var gedi_rh98_6k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-98-a0_vf_20190417_20230316_6000m")
//var gedi_rh98_12k = ee.Image('projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_rh-98-a0_vf_20190417_20230316_12000m')
var gedi_fhd_1k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_fhd-pai-1m-a0_vf_20190417_20230316_1000m")
var gedi_pai_1k = ee.Image("projects/ee-gedibio/assets/gedi/agg-stat-uploads/gediv002_pai-a0_vf_20190417_20230316_1000m")


// input args 
var als_img = coco_chm
var als_res = 1
var als_name = 'coco'
var als_year = 2019
var met_max_guess = 30
var alschm_agg_perc = 98

var bin_min = -99
var bin_max = 120
var bin_steps = 74

// RH98 -99 120 74 (3m bins)
// RH50 -99 120 147 (1.5m bins)
// PAI 0 12 49 (0.25 bins)
// FHD 0 4 41 (0.1 bins) used for GEDI 
//  - Suma: 0 0.5 21 (need to use diff bins since FHD calc is diff)

var gedi_img = gedi_rh98_1k
var metric = 'RH98'
var gedi_res = 1000
var gedi_fp_diam = 25

var not_urban_mask = nlcd_not_urban_mask
var not_water_mask = nlcd_not_water_mask

var dem = usned_dem



// ----------------------
// ----- PROCESSING -----
// ----------------------
var gedi_proj = gedi_img.projection()

// reformat the input GEDI image
gedi_img = gedi_img.addBands(gedi_img.select(7).multiply(ee.Number(Math.PI).multiply(ee.Number(gedi_fp_diam).divide(2).pow(2)))
                                               .divide(ee.Number(gedi_res).pow(2))
                                               .multiply(100))
                           .rename(['gedi_'+metric+'_mean', 
                                    'gedi_'+metric+'_meanbse', 
                                    'gedi_'+metric+'_med', 
                                    'gedi_'+metric+'_sd', 
                                    'gedi_'+metric+'_iqr', 
                                    'gedi_'+metric+'_p95', 
                                    'gedi_'+metric+'_shan', 
                                    'gedi_'+metric+'_shotcount',
                                    'gedi_'+metric+'_perccover'])
                           //.updateMask(gedi_img.select(7).gte(10)) 

// calculate slope from a DEM 
var slope = ee.Terrain.slope(dem).rename('slope')

// get the extent of ALS coverage
var als_reg = als_img.gte(0).reduceToVectors({
  scale: 100,  
  crs: gedi_proj,
  maxPixels: 1e13
})

print("Shot count histogram:", ui.Chart.image.histogram({
  image: gedi_img.select(7), 
  region: als_reg, 
  scale: gedi_res, 
  maxBuckets: 50, 
  minBucketWidth: 16}).setOptions({
    title: 'GEDI '+metric+' '+' (at '+gedi_res+'m resolution) shot count',
    titleX: 'GEDI Shot Count', 
    titleY: 'Frequency',
    fontName: 'arial',
    fontSize: 16,
    legend: {position: 'none'}
  }))
var gedi_img_sc_percs = ee.Dictionary(gedi_img.select(7).rename('sc').reduceRegion({
  reducer: ee.Reducer.percentile(ee.List.sequence(10,100,10)), 
  geometry: als_reg, 
  scale: gedi_res, 
  crs: gedi_proj}))
  
var gedi_img_sc_percs_ls = ee.List(gedi_img_sc_percs.values()).sort().reverse()

Map.addLayer(gedi_img.select(7), {min: 0, max: ee.Number(gedi_img_sc_percs_ls.get(0)).getInfo()}, 'GEDI '+metric+' shot count')
Map.addLayer(gedi_img.select(0), {min: 0, max: met_max_guess, palette: palettes.matplotlib.viridis[7]}, 'GEDI '+metric+' mean')
Map.addLayer(als_img, {min: 0, max: met_max_guess, palette: palettes.matplotlib.viridis[7]}, 'ALS '+metric+' (orig. res.)')


// combine masks to apply to ALS
// disturbance in the year of ALS acquisition or after
var forest_loss = ee.Image("UMD/hansen/global_forest_change_2022_v1_10").select('lossyear')
var not_disturb_mask = forest_loss.gte((ee.Number(als_year).subtract(2000))).selfMask().unmask().eq(0)
var valid_mask = not_disturb_mask.and(not_water_mask).and(not_urban_mask)

// resample ALS to match diameter of GEDI footprints
if (als_res < 20){
var alschm_agg_gedifp = als_img.updateMask(als_img.lte(100).and(als_img.gte(0)))
                              .reduceResolution({reducer: ee.Reducer.percentile({percentiles: [alschm_agg_perc]}), 
                                                 bestEffort: true})
                              .reproject({crs: gedi_proj, scale: gedi_fp_diam})
                              .rename('als_agg')
                              .updateMask(valid_mask)

  var als_agg_gedifp = alschm_agg_gedifp
} else {
  var als_agg_gedifp = als_img.reproject({crs: gedi_proj, scale: gedi_fp_diam})
                              .rename('als_agg')
                              .updateMask(valid_mask)
}
Map.addLayer(valid_mask.and(als_agg_gedifp.gte(0)), null, 'Valid Mask used on ALS', 0)

// Make a covering grid
var cov = als_reg.geometry().coveringGrid({proj: gedi_proj, 
                                           scale: gedi_res})
                            .map(function(g){
                              return g.set({'pc': als_reg.geometry().contains(g.geometry(), 5)})
                            }).filter(ee.Filter.eq('pc', true))

// Filter out grids which have more than 10% urban+water+disturbance area 
var cov_valid_rr = cov.map(function(g){
  var valid_red = ee.Image.pixelArea()
                    .updateMask(valid_mask)
                    .rename('valid')
                    .reduceRegion({
                      reducer: ee.Reducer.sum(), 
                      geometry: g.geometry(), 
                      scale: 30, 
                      crs: gedi_proj})
  var valid_perc = ee.Number(valid_red.get('valid')).divide(ee.Number(g.area(2))).multiply(100)
  return g.set(ee.Dictionary({'valid_perc': valid_perc}))
}).filter(ee.Filter.gte('valid_perc', 90))


// Define a function to compute multiple reducers on the ALS image
var als_rr_mult_stats = function(img, geom){
  var red = ee.Image(img).reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.percentile([25,50,75,95]),"",true) 
                              .combine(ee.Reducer.stdDev(), "", true)
                              .combine(ee.Reducer.fixedHistogram({min: bin_min, 
                                                                  max: bin_max, 
                                                                  steps: bin_steps, 
                                                                  cumulative: false}), "", true),
    geometry: ee.Geometry(geom), 
    scale: 25, 
    crs: gedi_proj})
  var names = ee.List(ee.Dictionary(red).keys())
  
  var nonzero_bins = ee.List(ee.Array(red.get('als_agg_histogram')).slice(1,1).project([0]).toList()).filter(ee.Filter.gt('item',0))
  var nonzero_sum = nonzero_bins.reduce(ee.Reducer.sum())
  var nonzero_ct = ee.Number(nonzero_bins.reduce(ee.Reducer.count()))
  var frac_bins = ee.List(nonzero_bins.map(function(i){
    var frac = ee.Number(i).divide(nonzero_sum)
    var frac_log = frac.log()
    return frac.multiply(frac_log)
  })) // TODO: verify there are at least 2 bins
  var shan = ee.Number(frac_bins.reduce(ee.Reducer.sum())).multiply(-1)
  var shan_nz = ee.Algorithms.If(nonzero_ct.gt(1), shan, null)
  
  return ee.Dictionary(['grid_res', gedi_res,
                        'als_agg_perc', alschm_agg_perc,
                        'als_agg_mean', red.get('als_agg_mean'),
                        'als_agg_med', red.get('als_agg_p50'),
                        'als_agg_sd', red.get('als_agg_stdDev'),
                        'als_agg_iqr', ee.Number(red.get('als_agg_p75')).subtract(ee.Number(red.get('als_agg_p25'))),
                        'als_agg_p95', red.get('als_agg_p95'),
                        'als_agg_shan', shan_nz,
                        'slope_mean', red.get('slope_mean'),
                        'slope_sd', red.get('slope_stdDev'),
                        'slope_p95', red.get('slope_p95')])
}


var cov_gedi_als_rr = cov_valid_rr.map(function(g){
  var gedi_red = ee.Image(gedi_img)
                   .reduceRegion({
                     reducer: ee.Reducer.first(), 
                     geometry: g.geometry(), 
                     scale: gedi_res, 
                     crs: gedi_proj})
  
  var als_red = als_rr_mult_stats(als_agg_gedifp.addBands(slope), g.geometry())

  return g.set(gedi_red)
          .set(als_red)
})

Map.addLayer(cov_gedi_als_rr, {'color': 'pink'}, 'Covering Grid Results (Res='+gedi_res+'m)', 0, 0.6)


var scatter_metric = function(fc, stat_label, y_als_metric_name, x_gedi_metric_name){
  var chart = ui.Chart.feature.byFeature(fc.limit(2000), x_gedi_metric_name, [y_als_metric_name])
  .setChartType('ScatterChart')
  .setOptions({title: 'GEDI '+metric+' '+stat_label+' at '+gedi_res+' m resolution',
               titleX: 'GEDI (m)', 
               titleY: 'ALS (m)',
               fontName: 'arial',
               fontSize: 16,
               pointSize: 4,
               dataOpacity: 0.3,
               legend: {position: 'right'},
               hAxis: {viewWindow: {min: 0, max: met_max_guess}},
               vAxis: {viewWindow: {min: 0, max: met_max_guess}},
               trendlines: {
                 0: {  // add a trend line to the 1st series
                 type: 'linear',  // or 'polynomial', 'exponential'
                 color: 'black',
                 lineWidth: 3,
                 opacity: 0.6,
                 visibleInLegend: true,
                showR2: true
                 }
               }
  })
  print(chart)
}

scatter_metric(cov_gedi_als_rr, 'mean', 'als_agg_mean', 'gedi_'+metric+'_mean')
scatter_metric(cov_gedi_als_rr, 'median', 'als_agg_med', 'gedi_'+metric+'_med')
scatter_metric(cov_gedi_als_rr, 'std. dev.', 'als_agg_sd', 'gedi_'+metric+'_sd')
scatter_metric(cov_gedi_als_rr, 'IQR', 'als_agg_iqr', 'gedi_'+metric+'_iqr')
scatter_metric(cov_gedi_als_rr, '95th Percentile', 'als_agg_p95', 'gedi_'+metric+'_p95')
scatter_metric(cov_gedi_als_rr, 'Shannons H', 'als_agg_shan', 'gedi_'+metric+'_shan')



// ----------------------
// ----- EXPORTS --------
// ----------------------
Export.table.toDrive({
  collection: cov_gedi_als_rr.filter(ee.Filter.gt(ee.String('gedi_')
                                                    .cat(ee.String(metric))
                                                    .cat(ee.String('_mean')),
                                                  -10)),
  description: 'gedi_'+metric+'_'+gedi_res+'m_val_als_'+alschm_agg_perc+'p_'+als_name, 
  fileNamePrefix: 'gedi_'+metric+'_'+gedi_res+'m_val_als_'+alschm_agg_perc+'p_'+als_name, 
  fileFormat: 'CSV'})