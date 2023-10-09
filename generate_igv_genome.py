"""
    Generates a IGV genome file in the new JSON format from a template (including updating gene file)
"""
import argparse
import json
import operator
import os
import subprocess
import tempfile
import urllib.request
import gzip
import shutil

# global variables
genome_json = {}

"""
sorting order for non integer chromosomes
"""
sorting_order = {
    'X': 100, 
    'Y': 101, 
    'MT': 102,
    "chr11_KI270721v1_random": 1000, 
    "chr14_GL000009v2_random": 1001, 
    "chr14_GL000225v1_random": 1002, 
    "chr14_KI270722v1_random": 1003, 
    "chr14_GL000194v1_random": 1004, 
    "chr14_KI270723v1_random": 1005, 
    "chr14_KI270724v1_random": 1006, 
    "chr14_KI270725v1_random": 1007, 
    "chr14_KI270726v1_random": 1008, 
    "chr15_KI270727v1_random": 1009, 
    "chr16_KI270728v1_random": 1010, 
    "chr17_GL000205v2_random": 1011, 
    "chr17_KI270729v1_random": 1012, 
    "chr17_KI270730v1_random": 1013, 
    "chr1_KI270706v1_random": 1014, 
    "chr1_KI270707v1_random": 1015, 
    "chr1_KI270708v1_random": 1016, 
    "chr1_KI270709v1_random": 1017, 
    "chr1_KI270710v1_random": 1018, 
    "chr1_KI270711v1_random": 1019, 
    "chr1_KI270712v1_random": 1020, 
    "chr1_KI270713v1_random": 1021, 
    "chr1_KI270714v1_random": 1022, 
    "chr22_KI270731v1_random": 1023, 
    "chr22_KI270732v1_random": 1024, 
    "chr22_KI270733v1_random": 1025, 
    "chr22_KI270734v1_random": 1026, 
    "chr22_KI270735v1_random": 1027, 
    "chr22_KI270736v1_random": 1028, 
    "chr22_KI270737v1_random": 1029, 
    "chr22_KI270738v1_random": 1030, 
    "chr22_KI270739v1_random": 1031, 
    "chr2_KI270715v1_random": 1032, 
    "chr2_KI270716v1_random": 1033, 
    "chr3_GL000221v1_random": 1034, 
    "chr4_GL000008v2_random": 1035, 
    "chr5_GL000208v1_random": 1036, 
    "chr9_KI270717v1_random": 1037, 
    "chr9_KI270718v1_random": 1038, 
    "chr9_KI270719v1_random": 1039, 
    "chr9_KI270720v1_random": 1040, 
    "chr1_KI270762v1_alt": 1041, 
    "chr1_KI270766v1_alt": 1042, 
    "chr1_KI270760v1_alt": 1043, 
    "chr1_KI270765v1_alt": 1044, 
    "chr1_GL383518v1_alt": 1045, 
    "chr1_GL383519v1_alt": 1046, 
    "chr1_GL383520v2_alt": 1047, 
    "chr1_KI270764v1_alt": 1048, 
    "chr1_KI270763v1_alt": 1049, 
    "chr1_KI270759v1_alt": 1050, 
    "chr1_KI270761v1_alt": 1051, 
    "chr2_KI270770v1_alt": 1052, 
    "chr2_KI270773v1_alt": 1053, 
    "chr2_KI270774v1_alt": 1054, 
    "chr2_KI270769v1_alt": 1055, 
    "chr2_GL383521v1_alt": 1056, 
    "chr2_KI270772v1_alt": 1057, 
    "chr2_KI270775v1_alt": 1058, 
    "chr2_KI270771v1_alt": 1059, 
    "chr2_KI270768v1_alt": 1060, 
    "chr2_GL582966v2_alt": 1061, 
    "chr2_GL383522v1_alt": 1062, 
    "chr2_KI270776v1_alt": 1063, 
    "chr2_KI270767v1_alt": 1064, 
    "chr3_JH636055v2_alt": 1065, 
    "chr3_KI270783v1_alt": 1066, 
    "chr3_KI270780v1_alt": 1067, 
    "chr3_GL383526v1_alt": 1068, 
    "chr3_KI270777v1_alt": 1069, 
    "chr3_KI270778v1_alt": 1070, 
    "chr3_KI270781v1_alt": 1071, 
    "chr3_KI270779v1_alt": 1072, 
    "chr3_KI270782v1_alt": 1073, 
    "chr3_KI270784v1_alt": 1074, 
    "chr4_KI270790v1_alt": 1075, 
    "chr4_GL383528v1_alt": 1076, 
    "chr4_KI270787v1_alt": 1077, 
    "chr4_GL000257v2_alt": 1078, 
    "chr4_KI270788v1_alt": 1079, 
    "chr4_GL383527v1_alt": 1080, 
    "chr4_KI270785v1_alt": 1081, 
    "chr4_KI270789v1_alt": 1082, 
    "chr4_KI270786v1_alt": 1083, 
    "chr5_KI270793v1_alt": 1084, 
    "chr5_KI270792v1_alt": 1085, 
    "chr5_KI270791v1_alt": 1086, 
    "chr5_GL383532v1_alt": 1087, 
    "chr5_GL949742v1_alt": 1088, 
    "chr5_KI270794v1_alt": 1089, 
    "chr5_GL339449v2_alt": 1090, 
    "chr5_GL383530v1_alt": 1091, 
    "chr5_KI270796v1_alt": 1092, 
    "chr5_GL383531v1_alt": 1093, 
    "chr5_KI270795v1_alt": 1094, 
    "chr6_GL000250v2_alt": 1095, 
    "chr6_KI270800v1_alt": 1096, 
    "chr6_KI270799v1_alt": 1097, 
    "chr6_GL383533v1_alt": 1098, 
    "chr6_KI270801v1_alt": 1099, 
    "chr6_KI270802v1_alt": 1100, 
    "chr6_KB021644v2_alt": 1101, 
    "chr6_KI270797v1_alt": 1102, 
    "chr6_KI270798v1_alt": 1103, 
    "chr7_KI270804v1_alt": 1104, 
    "chr7_KI270809v1_alt": 1105, 
    "chr7_KI270806v1_alt": 1106, 
    "chr7_GL383534v2_alt": 1107, 
    "chr7_KI270803v1_alt": 1108, 
    "chr7_KI270808v1_alt": 1109, 
    "chr7_KI270807v1_alt": 1110, 
    "chr7_KI270805v1_alt": 1111, 
    "chr8_KI270818v1_alt": 1112, 
    "chr8_KI270812v1_alt": 1113, 
    "chr8_KI270811v1_alt": 1114, 
    "chr8_KI270821v1_alt": 1115, 
    "chr8_KI270813v1_alt": 1116, 
    "chr8_KI270822v1_alt": 1117, 
    "chr8_KI270814v1_alt": 1118, 
    "chr8_KI270810v1_alt": 1119, 
    "chr8_KI270819v1_alt": 1120, 
    "chr8_KI270820v1_alt": 1121, 
    "chr8_KI270817v1_alt": 1122, 
    "chr8_KI270816v1_alt": 1123, 
    "chr8_KI270815v1_alt": 1124, 
    "chr9_GL383539v1_alt": 1125, 
    "chr9_GL383540v1_alt": 1126, 
    "chr9_GL383541v1_alt": 1127, 
    "chr9_GL383542v1_alt": 1128, 
    "chr9_KI270823v1_alt": 1129, 
    "chr10_GL383545v1_alt": 1130, 
    "chr10_KI270824v1_alt": 1131, 
    "chr10_GL383546v1_alt": 1132, 
    "chr10_KI270825v1_alt": 1133, 
    "chr11_KI270832v1_alt": 1134, 
    "chr11_KI270830v1_alt": 1135, 
    "chr11_KI270831v1_alt": 1136, 
    "chr11_KI270829v1_alt": 1137, 
    "chr11_GL383547v1_alt": 1138, 
    "chr11_JH159136v1_alt": 1139, 
    "chr11_JH159137v1_alt": 1140, 
    "chr11_KI270827v1_alt": 1141, 
    "chr11_KI270826v1_alt": 1142, 
    "chr12_GL877875v1_alt": 1143, 
    "chr12_GL877876v1_alt": 1144, 
    "chr12_KI270837v1_alt": 1145, 
    "chr12_GL383549v1_alt": 1146, 
    "chr12_KI270835v1_alt": 1147, 
    "chr12_GL383550v2_alt": 1148, 
    "chr12_GL383552v1_alt": 1149, 
    "chr12_GL383553v2_alt": 1150, 
    "chr12_KI270834v1_alt": 1151, 
    "chr12_GL383551v1_alt": 1152, 
    "chr12_KI270833v1_alt": 1153, 
    "chr12_KI270836v1_alt": 1154, 
    "chr13_KI270840v1_alt": 1155, 
    "chr13_KI270839v1_alt": 1156, 
    "chr13_KI270843v1_alt": 1157, 
    "chr13_KI270841v1_alt": 1158, 
    "chr13_KI270838v1_alt": 1159, 
    "chr13_KI270842v1_alt": 1160, 
    "chr14_KI270844v1_alt": 1161, 
    "chr14_KI270847v1_alt": 1162, 
    "chr14_KI270845v1_alt": 1163, 
    "chr14_KI270846v1_alt": 1164, 
    "chr15_KI270852v1_alt": 1165, 
    "chr15_KI270851v1_alt": 1166, 
    "chr15_KI270848v1_alt": 1167, 
    "chr15_GL383554v1_alt": 1168, 
    "chr15_KI270849v1_alt": 1169, 
    "chr15_GL383555v2_alt": 1170, 
    "chr15_KI270850v1_alt": 1171, 
    "chr16_KI270854v1_alt": 1172, 
    "chr16_KI270856v1_alt": 1173, 
    "chr16_KI270855v1_alt": 1174, 
    "chr16_KI270853v1_alt": 1175, 
    "chr16_GL383556v1_alt": 1176, 
    "chr16_GL383557v1_alt": 1177, 
    "chr17_GL383563v3_alt": 1178, 
    "chr17_KI270862v1_alt": 1179, 
    "chr17_KI270861v1_alt": 1180, 
    "chr17_KI270857v1_alt": 1181, 
    "chr17_JH159146v1_alt": 1182, 
    "chr17_JH159147v1_alt": 1183, 
    "chr17_GL383564v2_alt": 1184, 
    "chr17_GL000258v2_alt": 1185, 
    "chr17_GL383565v1_alt": 1186, 
    "chr17_KI270858v1_alt": 1187, 
    "chr17_KI270859v1_alt": 1188, 
    "chr17_GL383566v1_alt": 1189, 
    "chr17_KI270860v1_alt": 1190, 
    "chr18_KI270864v1_alt": 1191, 
    "chr18_GL383567v1_alt": 1192, 
    "chr18_GL383570v1_alt": 1193, 
    "chr18_GL383571v1_alt": 1194, 
    "chr18_GL383568v1_alt": 1195, 
    "chr18_GL383569v1_alt": 1196, 
    "chr18_GL383572v1_alt": 1197, 
    "chr18_KI270863v1_alt": 1198, 
    "chr19_KI270868v1_alt": 1199, 
    "chr19_KI270865v1_alt": 1200, 
    "chr19_GL383573v1_alt": 1201, 
    "chr19_GL383575v2_alt": 1202, 
    "chr19_GL383576v1_alt": 1203, 
    "chr19_GL383574v1_alt": 1204, 
    "chr19_KI270866v1_alt": 1205, 
    "chr19_KI270867v1_alt": 1206, 
    "chr19_GL949746v1_alt": 1207, 
    "chr20_GL383577v2_alt": 1208, 
    "chr20_KI270869v1_alt": 1209, 
    "chr20_KI270871v1_alt": 1210, 
    "chr20_KI270870v1_alt": 1211, 
    "chr21_GL383578v2_alt": 1212, 
    "chr21_KI270874v1_alt": 1213, 
    "chr21_KI270873v1_alt": 1214, 
    "chr21_GL383579v2_alt": 1215, 
    "chr21_GL383580v2_alt": 1216, 
    "chr21_GL383581v2_alt": 1217, 
    "chr21_KI270872v1_alt": 1218, 
    "chr22_KI270875v1_alt": 1219, 
    "chr22_KI270878v1_alt": 1220, 
    "chr22_KI270879v1_alt": 1221, 
    "chr22_KI270876v1_alt": 1222, 
    "chr22_KI270877v1_alt": 1223, 
    "chr22_GL383583v2_alt": 1224, 
    "chr22_GL383582v2_alt": 1225, 
    "chrX_KI270880v1_alt": 1226, 
    "chrX_KI270881v1_alt": 1227, 
    "chr19_KI270882v1_alt": 1228, 
    "chr19_KI270883v1_alt": 1229, 
    "chr19_KI270884v1_alt": 1230, 
    "chr19_KI270885v1_alt": 1231, 
    "chr19_KI270886v1_alt": 1232, 
    "chr19_KI270887v1_alt": 1233, 
    "chr19_KI270888v1_alt": 1234, 
    "chr19_KI270889v1_alt": 1235, 
    "chr19_KI270890v1_alt": 1236, 
    "chr19_KI270891v1_alt": 1237, 
    "chr1_KI270892v1_alt": 1238, 
    "chr2_KI270894v1_alt": 1239, 
    "chr2_KI270893v1_alt": 1240, 
    "chr3_KI270895v1_alt": 1241, 
    "chr4_KI270896v1_alt": 1242, 
    "chr5_KI270897v1_alt": 1243, 
    "chr5_KI270898v1_alt": 1244, 
    "chr6_GL000251v2_alt": 1245, 
    "chr7_KI270899v1_alt": 1246, 
    "chr8_KI270901v1_alt": 1247, 
    "chr8_KI270900v1_alt": 1248, 
    "chr11_KI270902v1_alt": 1249, 
    "chr11_KI270903v1_alt": 1250, 
    "chr12_KI270904v1_alt": 1251, 
    "chr15_KI270906v1_alt": 1252, 
    "chr15_KI270905v1_alt": 1253, 
    "chr17_KI270907v1_alt": 1254, 
    "chr17_KI270910v1_alt": 1255, 
    "chr17_KI270909v1_alt": 1256, 
    "chr17_JH159148v1_alt": 1257, 
    "chr17_KI270908v1_alt": 1258, 
    "chr18_KI270912v1_alt": 1259, 
    "chr18_KI270911v1_alt": 1260, 
    "chr19_GL949747v2_alt": 1261, 
    "chr22_KB663609v1_alt": 1262, 
    "chrX_KI270913v1_alt": 1263, 
    "chr19_KI270914v1_alt": 1264, 
    "chr19_KI270915v1_alt": 1265, 
    "chr19_KI270916v1_alt": 1266, 
    "chr19_KI270917v1_alt": 1267, 
    "chr19_KI270918v1_alt": 1268, 
    "chr19_KI270919v1_alt": 1269, 
    "chr19_KI270920v1_alt": 1270, 
    "chr19_KI270921v1_alt": 1271, 
    "chr19_KI270922v1_alt": 1272, 
    "chr19_KI270923v1_alt": 1273, 
    "chr3_KI270924v1_alt": 1274, 
    "chr4_KI270925v1_alt": 1275, 
    "chr6_GL000252v2_alt": 1276, 
    "chr8_KI270926v1_alt": 1277, 
    "chr11_KI270927v1_alt": 1278, 
    "chr19_GL949748v2_alt": 1279, 
    "chr22_KI270928v1_alt": 1280, 
    "chr19_KI270929v1_alt": 1281, 
    "chr19_KI270930v1_alt": 1282, 
    "chr19_KI270931v1_alt": 1283, 
    "chr19_KI270932v1_alt": 1284, 
    "chr19_KI270933v1_alt": 1285, 
    "chr19_GL000209v2_alt": 1286, 
    "chr3_KI270934v1_alt": 1287, 
    "chr6_GL000253v2_alt": 1288, 
    "chr19_GL949749v2_alt": 1289, 
    "chr3_KI270935v1_alt": 1290, 
    "chr6_GL000254v2_alt": 1291, 
    "chr19_GL949750v2_alt": 1292, 
    "chr3_KI270936v1_alt": 1293, 
    "chr6_GL000255v2_alt": 1294, 
    "chr19_GL949751v2_alt": 1295, 
    "chr3_KI270937v1_alt": 1296, 
    "chr6_GL000256v2_alt": 1297, 
    "chr19_GL949752v1_alt": 1298, 
    "chr6_KI270758v1_alt": 1299, 
    "chr19_GL949753v2_alt": 1300, 
    "chr19_KI270938v1_alt": 1301, 
    "chrUn_KI270302v1": 1302, 
    "chrUn_KI270304v1": 1303, 
    "chrUn_KI270303v1": 1304, 
    "chrUn_KI270305v1": 1305, 
    "chrUn_KI270322v1": 1306, 
    "chrUn_KI270320v1": 1307, 
    "chrUn_KI270310v1": 1308, 
    "chrUn_KI270316v1": 1309, 
    "chrUn_KI270315v1": 1310, 
    "chrUn_KI270312v1": 1311, 
    "chrUn_KI270311v1": 1312, 
    "chrUn_KI270317v1": 1313, 
    "chrUn_KI270412v1": 1314, 
    "chrUn_KI270411v1": 1315, 
    "chrUn_KI270414v1": 1316, 
    "chrUn_KI270419v1": 1317, 
    "chrUn_KI270418v1": 1318, 
    "chrUn_KI270420v1": 1319, 
    "chrUn_KI270424v1": 1320, 
    "chrUn_KI270417v1": 1321, 
    "chrUn_KI270422v1": 1322, 
    "chrUn_KI270423v1": 1323, 
    "chrUn_KI270425v1": 1324, 
    "chrUn_KI270429v1": 1325, 
    "chrUn_KI270442v1": 1326, 
    "chrUn_KI270466v1": 1327, 
    "chrUn_KI270465v1": 1328, 
    "chrUn_KI270467v1": 1329, 
    "chrUn_KI270435v1": 1330, 
    "chrUn_KI270438v1": 1331, 
    "chrUn_KI270468v1": 1332, 
    "chrUn_KI270510v1": 1333, 
    "chrUn_KI270509v1": 1334, 
    "chrUn_KI270518v1": 1335, 
    "chrUn_KI270508v1": 1336, 
    "chrUn_KI270516v1": 1337, 
    "chrUn_KI270512v1": 1338, 
    "chrUn_KI270519v1": 1339, 
    "chrUn_KI270522v1": 1340, 
    "chrUn_KI270511v1": 1341, 
    "chrUn_KI270515v1": 1342, 
    "chrUn_KI270507v1": 1343, 
    "chrUn_KI270517v1": 1344, 
    "chrUn_KI270529v1": 1345, 
    "chrUn_KI270528v1": 1346, 
    "chrUn_KI270530v1": 1347, 
    "chrUn_KI270539v1": 1348, 
    "chrUn_KI270538v1": 1349, 
    "chrUn_KI270544v1": 1350, 
    "chrUn_KI270548v1": 1351, 
    "chrUn_KI270583v1": 1352, 
    "chrUn_KI270587v1": 1353, 
    "chrUn_KI270580v1": 1354, 
    "chrUn_KI270581v1": 1355, 
    "chrUn_KI270579v1": 1356, 
    "chrUn_KI270589v1": 1357, 
    "chrUn_KI270590v1": 1358, 
    "chrUn_KI270584v1": 1359, 
    "chrUn_KI270582v1": 1360, 
    "chrUn_KI270588v1": 1361, 
    "chrUn_KI270593v1": 1362, 
    "chrUn_KI270591v1": 1363, 
    "chrUn_KI270330v1": 1364, 
    "chrUn_KI270329v1": 1365, 
    "chrUn_KI270334v1": 1366, 
    "chrUn_KI270333v1": 1367, 
    "chrUn_KI270335v1": 1368, 
    "chrUn_KI270338v1": 1369, 
    "chrUn_KI270340v1": 1370, 
    "chrUn_KI270336v1": 1371, 
    "chrUn_KI270337v1": 1372, 
    "chrUn_KI270363v1": 1373, 
    "chrUn_KI270364v1": 1374, 
    "chrUn_KI270362v1": 1375, 
    "chrUn_KI270366v1": 1376, 
    "chrUn_KI270378v1": 1377, 
    "chrUn_KI270379v1": 1378, 
    "chrUn_KI270389v1": 1379, 
    "chrUn_KI270390v1": 1380, 
    "chrUn_KI270387v1": 1381, 
    "chrUn_KI270395v1": 1382, 
    "chrUn_KI270396v1": 1383, 
    "chrUn_KI270388v1": 1384, 
    "chrUn_KI270394v1": 1385, 
    "chrUn_KI270386v1": 1386, 
    "chrUn_KI270391v1": 1387, 
    "chrUn_KI270383v1": 1388, 
    "chrUn_KI270393v1": 1389, 
    "chrUn_KI270384v1": 1390, 
    "chrUn_KI270392v1": 1391, 
    "chrUn_KI270381v1": 1392, 
    "chrUn_KI270385v1": 1393, 
    "chrUn_KI270382v1": 1394, 
    "chrUn_KI270376v1": 1395, 
    "chrUn_KI270374v1": 1396, 
    "chrUn_KI270372v1": 1397, 
    "chrUn_KI270373v1": 1398, 
    "chrUn_KI270375v1": 1399, 
    "chrUn_KI270371v1": 1400, 
    "chrUn_KI270448v1": 1401, 
    "chrUn_KI270521v1": 1402, 
    "chrUn_GL000195v1": 1403, 
    "chrUn_GL000219v1": 1404, 
    "chrUn_GL000220v1": 1405, 
    "chrUn_GL000224v1": 1406, 
    "chrUn_KI270741v1": 1407, 
    "chrUn_GL000226v1": 1408, 
    "chrUn_GL000213v1": 1409, 
    "chrUn_KI270743v1": 1410, 
    "chrUn_KI270744v1": 1411, 
    "chrUn_KI270745v1": 1412, 
    "chrUn_KI270746v1": 1413, 
    "chrUn_KI270747v1": 1414, 
    "chrUn_KI270748v1": 1415, 
    "chrUn_KI270749v1": 1416, 
    "chrUn_KI270750v1": 1417, 
    "chrUn_KI270751v1": 1418, 
    "chrUn_KI270752v1": 1419, 
    "chrUn_KI270753v1": 1420, 
    "chrUn_KI270754v1": 1421, 
    "chrUn_KI270755v1": 1422, 
    "chrUn_KI270756v1": 1423, 
    "chrUn_KI270757v1": 1424, 
    "chrUn_GL000214v1": 1425, 
    "chrUn_KI270742v1": 1426, 
    "chrUn_GL000216v2": 1427, 
    "chrUn_GL000218v1": 1428, 
    "chrY_KI270740v1_random": 1429
    }

def parse_args():
    """
                parses the arguments

    :return: argparse object containing all provided arguments
    """

    print("parsing args...")

    parser = argparse.ArgumentParser(description="Generates a IGV genome JSON")
    parser.add_argument("template_file", help="Template JSON containing all links to the input files.")
    parser.add_argument("hgnc_file", help="file path to the the HGNC table file (containing HGNC id <-> gene name mapping")
    parser.add_argument("output", help="file path for the generated output IGV genome JSON file")

    return parser.parse_args()


def parse_json(template_file: str):
    """
            parses the template file

    :param template_file:   file path to the json template file
    """
    global genome_json
    with open(template_file, 'r') as json_file:
        genome_json = json.load(json_file)
    return


def download_files(output_folder: str):
    """
            downloads all distant files in the template json and links to them

    :param output_folder:   target folder for the downloads

    :return:             modified JSON template with links to the local files
    """
    global genome_json
    # download all required files

    for key in ["fastaURL", "indexURL", "cytobandURL", "aliasURL"]:
        url = genome_json[key]
        filename = os.path.basename(url)
        print("Downloading file '" + filename + "'...")
        urllib.request.urlretrieve(url, os.path.join(output_folder, filename))
        genome_json[key] = filename

    # download tracks
    for track in genome_json["tracks"]:
        url = track["url"]
        filename = os.path.basename(url)
        print("Downloading file '" + filename + "'...")
        urllib.request.urlretrieve(url, os.path.join(output_folder, filename))
        track["url"] = filename
        # load optional index
        if "indexURL" in track:
            url = track["indexURL"]
            filename = os.path.basename(url)
            print("Downloading file '" + filename + "'...")
            urllib.request.urlretrieve(url, os.path.join(output_folder, filename))
            track["indexURL"] = filename

    return


def load_hgnc_file(hgnc_filepath):
    print("Parsing HGNC file...")
    hgnc_mapping = {}
    with open(hgnc_filepath, 'r', encoding="utf8") as hgnc_file:
        for line in hgnc_file:
            if line.startswith("hgnc_id\tsymbol\tname"):
                # skip header
                continue
            elif line.startswith("HGNC:"):
                # parse line
                split_line = line.split('\t')
                hgnc_id = int(split_line[0].split(':')[1])
                symbol = split_line[1].strip()
                hgnc_mapping[hgnc_id] = symbol
    print("\t " + str(len(hgnc_mapping.keys())) + " ids parsed.")
    return hgnc_mapping


def sort_gff3(content):
    """
                sorts the content gff3 file based on the sorting order defined above

    :param content:  content of the gff file (without headers)
    :return:         sorted list of lists
    """

    print("sorting gff3 data...")

    # preprocess chromosome column and positioning column
    for line in content:
        if line[0] in sorting_order:
            line[0] = sorting_order[line[0]]
        elif line[0].startswith("GL0"):
            line[0] = 10000 + int(float(line[0][2:]) * 10) #move GL0000XX.X to the end
        elif line[0].startswith("KI"):
            line[0] = 20000 + int(float(line[0][2:]) * 10) #move KI0000XX.X to the end
        else:
            line[0] = int(line[0])

        line[3] = int(line[3])
        line[4] = int(line[4])

    # sort content
    content.sort(key=operator.itemgetter(0, 3))

    # create reverse sorting_order:
    reverse_sorting_order = {}
    for k, v in sorting_order.items():
        reverse_sorting_order[v] = k

    # replace placeholder with original chromosome names:
    for line in content:
        if line[0] in reverse_sorting_order:
            line[0] = reverse_sorting_order[line[0]]

    return content

def update_gene_file(hgnc_mapping, output_folder):
    """

    :param hgnc_mapping:
    :return:
    """

    print("Modifying GFF3 files (updating gene names) ...")
    global genome_json
    # modify all gff3 track files:
    for track in genome_json["tracks"]:
        if "format" in track and track["format"] == "gff3":

            print("Modifying GFF3 file '" + track["url"] + "'...")

            # unzip and modify
            n_comment_lines = 0
            n_unmodified_lines = 0
            n_modified_lines = 0
            n_ignored = 0
            comment_lines = []
            content_lines = []
            with gzip.open(os.path.join(output_folder, track["url"]), 'rb') as compressed_gff3:

                    for line in compressed_gff3:
                        line = line.decode('utf-8')
                        # skip comments
                        if line.startswith("#"):
                            if line.strip() == "###":
                                # ignore
                                n_ignored += 1
                                continue
                            comment_lines.append(line)
                            n_comment_lines += 1
                            continue
                        # detect HGNC ids
                        annotation_column = line.split('\t')[8]
                        if "[Source:HGNC Symbol%3BAcc:" in annotation_column:
                            kv_list = annotation_column.split(';')
                            idx_name = -1
                            idx_description = -1
                            for idx in range(len(kv_list)):
                                if kv_list[idx].startswith("Name="):
                                    idx_name = idx
                                elif kv_list[idx].startswith("description="):
                                    idx_description = idx
                            if idx_description > -1:
                                # extract HGNC id
                                hgnc_id = int(kv_list[idx_description].split('[')[1].split(']')[0].split(':')[-1])

                                if hgnc_id not in hgnc_mapping.keys():
                                    print("Warning: HGNC id " + str(hgnc_id) + " not found in HGNC file!")
                                    # store line unmodified
                                    content_lines.append(line.split('\t'))
                                    n_unmodified_lines += 1
                                    continue
                                if idx_name > -1:
                                    kv_list[idx_name] = "Name=" + hgnc_mapping[hgnc_id]
                                else:
                                    kv_list.append("Name=" + hgnc_mapping[hgnc_id])

                                split_line = line.split('\t')
                                split_line[8] = ";".join(kv_list)
                                content_lines.append(split_line)
                                n_modified_lines += 1
                        else:
                            # no HGNC identifier -> store unmodified
                            content_lines.append(line.split('\t'))
                            n_unmodified_lines += 1

            # stats
            print("\tcomment lines: " + str(n_comment_lines))
            print("\tunmodified lines: " + str(n_unmodified_lines))
            print("\tmodified lines: " + str(n_modified_lines))
            print("\tignored lines: " + str(n_ignored))

            # sort content

            content_lines = sort_gff3(content_lines)

            # write modified file to disk
            print("Writing modified file to disk...")
            with open(os.path.join(output_folder, os.path.splitext(track["url"])[0]), 'wt') as modified_gff3:
                for line in comment_lines:
                    modified_gff3.write(line)
                for line in content_lines:
                    modified_gff3.write("\t".join(str(e) for e in line))

            # bgzip
            print("Compressing file...")
            print("\tignored lines: " + str(n_ignored))
            rc = subprocess.call(["bgzip", "-f", os.path.join(output_folder, os.path.splitext(track["url"])[0])])
            if rc != 0:
                raise RuntimeError("bgzip failed with return code " + str(rc) + "!")

            # tabix
            print("Indexing file...")
            rc = subprocess.call(["tabix", "-p", "gff", os.path.join(output_folder, track["url"])])
            if rc != 0:
                raise RuntimeError("tabix failed with return code " + str(rc) + "!")
            # add index to JSON
            track["indexURL"] = track["url"] + ".tbi"
    return


def update_alias_file(output_folder):
    """
                Extends the alias tab file
    """
    print("Updating alias file...")
    alias_file_name = genome_json["aliasURL"]
    file_buffer = []
    with open(os.path.join(output_folder, alias_file_name), 'r') as alias_file:
        for line in alias_file:
            if line.startswith("chrM"):
                # add 'chrMT' and 'M' as valid aliases
                line = line.strip() + "\tchrMT\tM\n"
            file_buffer.append(line)

    with open(os.path.join(output_folder, alias_file_name), 'w') as alias_file:
        alias_file.writelines(file_buffer)
    return


def main():

    args = parse_args()

    # read template
    parse_json(args.template_file)

    # download files to local storage
    output_folder = os.path.dirname(args.output)
    download_files(output_folder)

    # load hgnc file
    hgnc_mapping = load_hgnc_file(args.hgnc_file)

    # update gene file
    update_gene_file(hgnc_mapping, output_folder)

    # update alias file
    update_alias_file(output_folder)

    # store modified JSON file
    with open(args.output, 'w') as output_file:
        print("Writing genome JSON file...")
        json.dump(genome_json, output_file, indent=4)

    print("\nfinished.")


if __name__ == '__main__':
    main()
