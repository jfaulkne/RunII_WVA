# $Id: runMakeTreeFromPAT_cfg.py,v 1.9 2013/01/24 15:42:54 mschrode Exp $
#
# Expects a file name as argument e.g.
# cmsRun runMakeLostLeptonTreeFromPAT_cfg.py dataset=/store/user/mschrode/HT/RA2PreSelection_Run2012A-13Jul2012-v1_V4/21a074f94cdbe7cfbeeb19be46b40a6a/RA2Skim_9_1_h6A.root
# cmsRun ../test/runMakeLostLeptonTreeFromPAT_cfg.py dataset=/store/user/mschrode/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/RA2PreSelection_Summer12_DR53X-PU_S10_START53_V7A-v1_V4/6c50609e978ba7d5388d5439fc628605/RA2Skim_100_1_dgv.root, global_tag=START53_V7F::All, MC=True, debug=True

# Read parameters
from SemiLeptonic.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()

MC= parameters.value("MC", False)
CH = parameters.value("channel","")
RUN = parameters.value("run",0)
ERA = parameters.value("era","")
OUTFILE = parameters.value("outfile","ReducedSelection")

SingleLepton_Dataset = {
"el":{"prompt":[
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/604/00000/AE22AF42-902A-E511-8A22-02163E012B30.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/612/00000/50DA7894-932A-E511-801E-02163E0136A2.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/628/00000/40EF63A0-B52A-E511-8B57-02163E0133F0.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/638/00000/0CDDB666-E72A-E511-9BFD-02163E011DE5.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/638/00000/B2FC1038-372B-E511-AA94-02163E013481.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/642/00000/2C622272-D02A-E511-9F20-02163E013645.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/0E4B7E28-8D2C-E511-BFDA-02163E01477B.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/1247CF12-932C-E511-B9ED-02163E01354D.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/3A437BCB-912C-E511-96D0-02163E012934.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/7077210E-8F2C-E511-97D5-02163E0138EC.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/9EFCB7EB-C12C-E511-B8BB-02163E012158.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/A42F5B12-9C2C-E511-AAB3-02163E0134C3.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/C2E62796-942C-E511-8869-02163E0121A1.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/CCA6600A-912C-E511-B1EF-02163E0133D0.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/D84DA8FC-8F2C-E511-9B5A-02163E01420D.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/716/00000/02C46BE4-302C-E511-A01C-02163E0128BF.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/721/00000/9663EE89-022C-E511-985A-02163E01359E.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/721/00000/CADC920F-E02B-E511-BB9B-02163E01412F.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/781/00000/662647FF-9B2C-E511-85C5-02163E012965.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/00CD59FD-2B2D-E511-8DB2-02163E01267F.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/500D754A-292D-E511-AEBA-02163E0118F2.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/CA2F0F5F-242D-E511-96AB-02163E011D88.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/CAA3B776-312D-E511-923A-02163E01280D.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/D030EA86-9C2D-E511-8FCB-02163E014181.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/102/00000/06E11C1F-DA2F-E511-B04A-02163E0117FF.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/116/00000/2EE8BBD0-7730-E511-95F5-02163E013560.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/126/00000/7E4ACB62-6E31-E511-858C-02163E0133F0.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/126/00000/AA0C0594-CA30-E511-BD15-02163E0123C0.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/488/00000/B4D62F1C-EE35-E511-BAFE-02163E01477B.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/489/00000/686DE8AC-3334-E511-BAF8-02163E0144F6.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/490/00000/E2682AAA-3334-E511-9706-02163E012B10.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/492/00000/16665670-4734-E511-9940-02163E013901.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/496/00000/2CA28248-CC35-E511-BB41-02163E01258B.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/499/00000/DC5562FB-CB35-E511-99EB-02163E0125D6.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/501/00000/A27F94C6-CA35-E511-82EC-02163E01459D.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/571/00000/1C77574A-F13D-E511-9CD6-02163E013473.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/578/00000/78DD4869-F63D-E511-B425-02163E01289E.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/596/00000/821E01AC-053E-E511-8B0B-02163E01289E.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/607/00000/863448CE-0E3E-E511-B24F-02163E0128F0.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/608/00000/FE15468F-0F3E-E511-9A77-02163E011B3C.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/611/00000/B8127DE6-103E-E511-A5FB-02163E013542.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/612/00000/C0E78833-123E-E511-9E55-02163E014316.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/253/620/00000/C21CA147-763F-E511-A201-02163E0143B7.root'
],"july":[
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/109023E8-C02E-E511-8913-0025905A610A.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/18EA30CE-C02E-E511-AAAE-0025905B859E.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/1C847DE8-C02E-E511-8D41-0025905A60DA.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/3861A09F-C02E-E511-B4B3-0025905A613C.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/3A979A2B-C12E-E511-A5F7-0025905A48D6.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/3C61DCE6-C02E-E511-ADD1-002618943971.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/3E4B2930-C12E-E511-872B-002590593920.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/48C83128-C12E-E511-B235-003048FFCC0A.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/4A52C8E2-C02E-E511-8FF3-0025905A60A0.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/5097CCE8-C02E-E511-B31C-002618FDA28E.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/52D4A025-C12E-E511-89FD-00261894382D.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/584065EA-C02E-E511-A27B-0025905A613C.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/58C331BE-C02E-E511-A2B2-0025905A608C.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/6012E4D2-C02E-E511-A24C-0025905B858E.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/6CFFA52B-C12E-E511-B60B-0025905A60F2.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/783F7796-C02E-E511-ABBE-003048FFD770.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/7E1814E8-C02E-E511-8E7E-0025905AA9CC.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/9261D028-C12E-E511-BBBF-0025905A60B4.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/96D4049C-C02E-E511-92C9-00261894397A.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/B2000C33-C12E-E511-BD6C-002618943849.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/B6D0E4CB-C12E-E511-AA9F-0025905B85AA.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/C68F6A2C-C12E-E511-8407-0025905A48BC.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/CC503EEC-C02E-E511-93A3-002618943982.root',
'/store/data/Run2015B/SingleElectron/MINIAOD/17Jul2015-v1/30000/F40EF063-C02E-E511-A074-0025905A612A.root'
],"254833":[
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/027C5D36-524B-E511-B191-02163E014705.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/04D13977-334B-E511-BD56-02163E0121EB.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/1C60B365-334B-E511-9563-02163E014742.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/3A56A861-334B-E511-B711-02163E011D3C.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/561E0D69-334B-E511-AB10-02163E0139DB.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/66446468-334B-E511-8642-02163E013440.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/7C8B3F67-334B-E511-8D28-02163E0137B2.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/8E23446D-334B-E511-90DF-02163E0126FB.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/9A46897B-334B-E511-A570-02163E014615.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/AA78446F-334B-E511-B5D8-02163E013818.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/C82CED39-524B-E511-8790-02163E014705.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/DCDF5B3B-524B-E511-B9B2-02163E014705.root',
'/store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/833/00000/EA0A1F7C-334B-E511-81D8-02163E01227F.root'
]},
"mu":{"prompt":[
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/604/00000/1606A3BF-972A-E511-86A7-02163E013515.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/612/00000/7A0CE8FF-A72A-E511-B7DC-02163E011D1C.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/638/00000/10C07DF0-FA2A-E511-846A-02163E01474A.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/642/00000/D20F8A8A-DE2A-E511-9D16-02163E0133FF.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/643/00000/1E72D617-BE2C-E511-96A0-02163E0139A2.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/643/00000/3C563818-BE2C-E511-993B-02163E0144D6.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/643/00000/C061E81E-BE2C-E511-AA5F-02163E01208E.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/643/00000/C6F9761A-BE2C-E511-932C-02163E011D30.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/643/00000/CC199B16-BE2C-E511-B1A5-02163E012B30.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/643/00000/D4FE721A-BE2C-E511-856B-02163E01250C.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/721/00000/30BD3F39-542C-E511-A5DA-02163E011DA4.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/721/00000/329A9A3B-542C-E511-B6B8-02163E01360E.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/781/00000/CE3105FD-A82C-E511-B330-02163E011955.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/883/00000/089D049E-262D-E511-85A7-02163E0146EB.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/883/00000/62919ECB-1F2D-E511-B387-02163E013796.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/883/00000/E2546D9E-492D-E511-9977-02163E011D46.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/883/00000/E49261AB-492D-E511-9FCA-02163E011E24.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/116/00000/7E03BD9B-7730-E511-BA13-02163E011976.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/201D1572-EB31-E511-852A-02163E01340A.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/2471A7F2-6C31-E511-A71B-02163E0133F0.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/58D01EE9-6B31-E511-852A-02163E0139B0.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/647C27E8-E231-E511-A3D9-02163E012704.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/9E0526D6-F131-E511-94B7-02163E01428F.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/A0BF6FF9-8231-E511-AD78-02163E0133AD.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/BA9B25ED-E231-E511-976A-02163E01340A.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/BE4EB725-7131-E511-84DD-02163E0139B0.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/126/00000/DAF3D992-ED31-E511-812B-02163E014181.root'
], "july":[
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/16B50792-172E-E511-B0C8-0025905C43EC.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/18DDE791-172E-E511-82A7-0025905C3D6A.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/20A12A92-172E-E511-9EBB-0025904C68DC.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/2A679E94-172E-E511-B36C-0025905C96A4.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/32CA1C92-172E-E511-BD1F-0025904C68DC.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/3A314794-172E-E511-BACB-0025905C2CB8.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/5AC60C92-172E-E511-8161-0025904CF93C.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/685CF492-172E-E511-82C3-0025905C96A4.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/760EDD94-172E-E511-B801-0025905C2D9A.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/76384094-172E-E511-833C-0025905C2D98.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/7E1D2E96-172E-E511-B304-0025904CF100.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/7ED0F296-172E-E511-A7D8-0025905C94D2.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/8835A890-172E-E511-9456-0025905C2CBA.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/8C8B7A91-172E-E511-B582-0025905C3DF6.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/960C8C77-172E-E511-8303-0025904C66E6.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/96C17592-172E-E511-98F9-0025904B7C40.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/A8016493-172E-E511-8A21-0025905C42B8.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/AC31B090-172E-E511-A9CF-0025905C2CE8.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/BEBDE260-172E-E511-97AA-0025905C3D3E.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/D0C57E91-172E-E511-8B70-0025905BA736.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/D43EF291-172E-E511-99D2-0025905C3D6A.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/E43DA791-172E-E511-BB8B-0025905BA736.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/FC24B692-172E-E511-8421-0025905C2CBE.root',
'/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/FC32DF92-172E-E511-AD49-0025905C42B6.root'
],"254833":[
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/220E01C3-104B-E511-837F-02163E015541.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/4055C6BB-104B-E511-B242-02163E0121C9.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/640CD4BE-104B-E511-B7CE-02163E011EF4.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/68B440BA-104B-E511-80B5-02163E0133A9.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/8AEF9BBB-104B-E511-8EF2-02163E0124FE.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/A20B94B7-104B-E511-8341-02163E013964.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/BE9F48BC-104B-E511-ADDD-02163E0141ED.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/EAEF3DB9-104B-E511-9C21-02163E0139BA.root',
'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/FCB9CBAE-104B-E511-BF89-02163E011F65.root'
]}}

if MC: dataSetName = parameters.value("dataset","")
else: dataSetName = parameters.value("dataset", SingleLepton_Dataset[CH][ERA][RUN])

JSON = {"prompt":"/afs/cern.ch/work/j/jfaulkne/CMSSW_7_4_7_patch2/src/SemiLeptonicWVA/json/Cert_246908-254879_13TeV_PromptReco_Collisions15_JSON_prompt.txt",
"july":"/afs/cern.ch/work/j/jfaulkne/CMSSW_7_4_7_patch2/src/SemiLeptonicWVA/json/Cert_246908-254879_13TeV_PromptReco_Collisions15_JSON_july.txt",
"254833":"/afs/cern.ch/work/j/jfaulkne/CMSSW_7_4_7_patch2/src/SemiLeptonicWVA/json/Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt"}

global_tag = parameters.value("global_tag","74X_dataRun2_v2")
QCD= parameters.value("QCD", False)
LostLepton= parameters.value("LostLepton", False)
debug= parameters.value("debug", False)
nJetsMin    = parameters.value("njets_min",0)
htMin       = parameters.value("ht_min",0)
mhtMin      = parameters.value("mht_min",0)
NumProcessedEvt=parameters.value("NumProcessedEvt",-1)
METFiltersProcess=parameters.value("METFiltersProcess","")
DoAK8Reclustering=parameters.value("DoAK8Reclustering",True)
DoJECCorrection=parameters.value("DoJECCorrection",True)
DoPuppi=parameters.value("DoPuppi",False)
LeptonFilter=parameters.value("leptonFilter",True)
GenJetsAK8Reclustering=parameters.value("genJetsAK8Reclustering",True)
isHBHEEarlyData = parameters.value("isHBHEEarlyData",True)
JsonFileName=parameters.value("jsonFileName",JSON[ERA])
IsCrab=parameters.value("isCrab",False)

processName      = parameters.value("name","RSGraviton1000")


print "***** SETUP ************************************"
print "  dataSetName : "+dataSetName
print " global_tag : "+global_tag
print " runningOnMC : "+str(MC)
print " runningOnQCD : "+str(QCD)
print " LostLepton(MC) : "+str(LostLepton)
print "     nJetsMin : "+str(nJetsMin)
print "        htMin : "+str(htMin)
print "       mhtMin : "+str(mhtMin)
print "       debug : "+str(debug)
print "       num of events : "+str(NumProcessedEvt)
print "       doAK8Reclustering : "+str(DoAK8Reclustering)
print "       doJECCorrection : "+str(DoJECCorrection)
print "       doPuppi : "+str(DoPuppi)
print "       leptonFilter : "+str(LeptonFilter)
print "       genJetsAK8Reclustering : "+str(GenJetsAK8Reclustering)
print "       isHBHEEarlyData : "+str(isHBHEEarlyData)
print "       jsonFileName : "+str(JsonFileName)
print "       isCrab : "+str(False)
print "************************************************"

# The process needs to be defined AFTER reading sys.argv,
# otherwise edmConfigHash fails
import FWCore.ParameterSet.Config as cms
#process = cms.Process("RA2EventSelection")
process = cms.Process(processName)

from SemiLeptonicWVA.TreeMaker.makeTreeFromMiniAOD_cff import makeTreeTreeFromMiniAOD
makeTreeTreeFromMiniAOD(process,
                outFileName=OUTFILE,
                NJetsMin=nJetsMin,
                HTMin=htMin,
                MHTMin=mhtMin,
                reportEveryEvt=1000,
                testFileName=dataSetName,
		Global_Tag=global_tag,
                METFiltersProcess=METFiltersProcess,
		MC=MC,
		QCD=QCD,
		LostLepton=LostLepton,
		debug = debug,
                numProcessedEvt=NumProcessedEvt,
                doAK8Reclustering=DoAK8Reclustering,
                doJECCorrection=DoJECCorrection,
                doPuppi=DoPuppi,
                leptonFilter=LeptonFilter,
                genJetsAK8Reclustering=GenJetsAK8Reclustering,
                customizeHBHENoiseForEarlyData=isHBHEEarlyData,
                jsonFileName=JsonFileName,
                isCrab=IsCrab)

