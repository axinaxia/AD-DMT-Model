library(RMySQL)
con <- dbConnect(MySQL(), 
                 user='linus.jonsson', 
                 password=rstudioapi::askForPassword(prompt = "Password"),
                 host='h1cogbase01.nvs.ki.se', 
                 port=3306,
                 dbname='priv_linjon')
#dbExecute(con, 'SET GLOBAL local_infile = "ON"')
#Create tables on server----

dbExecute(con, 'DROP TABLE IF EXISTS priv_linjon.INCLUSION;')
dbExecute(con, 
          "CREATE TABLE priv_linjon.INCLUSION AS 
  SELECT t.*, b.amyloid FROM (
    SELECT LOPNR, DIAGNOS, d_YOB, DIAGNOSDATUM, ROW_NUMBER() OVER (PARTITION BY LOPNR ORDER BY DIAGNOSDATUM) AS n
    FROM SVEDEM_2021.SVEDEM s
    WHERE SOURCE_TABLE IN ('GRUNDREG__2', 'GRUNDREG__3')
  ) t 
      LEFT JOIN priv_linjon.BM_PATH b on t.LOPNR=b.LOPNR
  WHERE n=1")

dbExecute(con, 'DROP TABLE IF EXISTS priv_linjon.DEMO;')
dbExecute(con, 
          "CREATE TABLE priv_linjon.DEMO AS 
  SELECT * FROM (
    SELECT *, ROW_NUMBER() OVER (PARTITION BY LOPNR ORDER BY DIAGNOSDATUM) AS n
    FROM SVEDEM_2021.SVEDEM s
    WHERE SOURCE_TABLE IN ('GRUNDREG__2', 'GRUNDREG__3')
  ) t 
  WHERE n=1")

dbExecute(con, 'DROP TABLE IF EXISTS priv_linjon.COG;')

dbExecute(con, 
          "CREATE TABLE priv_linjon.COG AS 
  SELECT LOPNR, DIAGNOSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.GRUNDREG__2
          UNION ALL
          SELECT LOPNR, DIAGNOSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.GRUNDREG__3
          UNION ALL
          SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.UPPF__2
          UNION ALL
          SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, MOCA, MOCA_VARDE
          FROM SVEDEM_2021.UPPF__3
          UNION ALL
          SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, NULL AS MOCA, NULL AS MOCA_VARDE
          FROM SVEDEM_2021.UPPF_HE_2
          UNION ALL
          SELECT LOPNR, UPPFOLJNINGSDATUM as DATE, MMSESR, MMSESR_VARDE, NULL AS MOCA, NULL AS MOCA_VARDE
          FROM SVEDEM_2021.UPPF_SA_2
          ;")

dbExecute(con, 'DROP TABLE IF EXISTS priv_linjon.MORT;')
dbExecute(con, 
          "CREATE TABLE priv_linjon.MORT AS 
          
            SELECT LOPNR, DODSDAT
            FROM SVEDEM_2021.DORS
            UNION ALL
            SELECT LOPNR, DODSDAT
            FROM SVEDEM_2021.DORS_AVI
          ;")

dbExecute(con, 'DROP TABLE IF EXISTS priv_linjon.INST;')
dbExecute(con, 
          "CREATE TABLE priv_linjon.INST AS 
            SELECT LOPNR, min(PERIOD) AS PERIOD 
            FROM SVEDEM_2021.SOL
            WHERE BOFORM=2
            GROUP BY LOPNR
          ;")

# Vitamin K antagonists (ATC B01AA)
# Heparin, LMWH and related agents (ATC B01AB, B01AX05)
# Other antiplatelet agents (ATC B01AC)
# Direct oral anticoagulants (ATC B01AE, B01AF)

DROP TABLE IF EXISTS priv_linjon.PHARM;
CREATE TABLE priv_linjon.PHARM AS 
SELECT S.LOPNR, ATC, subnamn, EDATUM
FROM SVEDEM_2021.LMED S
WHERE (substr(ATC,1,5) IN ('B01AA', 'B01AB', 'B01AC', 'B01AE', 'B01AF') OR substr(ATC,1,7)='B01AX05');

CREATE TABLE priv_linjon.ANTICOAG AS 
SELECT p.LOPNR, ATC, subnamn, min(EDATUM) as first_date, max(EDATUM) as last_date
FROM priv_linjon.PHARM p
JOIN priv_linjon.INCLUSION i on p.LOPNR=i.LOPNR 
GROUP BY p.LOPNR, ATC, subnamn;
