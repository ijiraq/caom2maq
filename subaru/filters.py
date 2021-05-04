from astropy import units

DISPERSERS = {"75": {"SY47": (4700, 9100, 250),
                     "SO58": (5800, 10000, 250),
                     "DEFAULT": (3000, 11000, 250)},
              "150": {"L550": (3400, 5500, 700),
                      "SO58": (5800, 10000, 500),
                      "SY47": (4700, 9100, 500),
                      "DEFAULT": (3000, 11000, 600)},
              "300B": {"SY47": (4700, 9100, 1000),
                       "L600": (3700, 6000, 1000),
                       "DEFAULT": (3000, 11000, 1000)},
              "300R": {"SY47": (4900, 9100, 1000),
                       "SO58": (5800, 10000, 1000),
                       "L600": (3700, 5950, 2000),
                       "L550": (3400, 5250, 2000),
                       "DEFAULT": (3000, 11000, 1500)},
              "Echelle": {"I": (7300, 8700, 2500),
                          "z": (8300, 10000, 2500),
                          "DEFAULT": (3000, 11000, 2500)},
              "VPH450": {"NONE": (3800, 5250, 3000),
                         "DEFAULT":(3800, 5250, 3000)},
              "VPH520": {"NONE": (4450, 6050, 3000),
                         "DEFAULT": (4450, 6050, 3000)},
              "VPH650": {"SY47": (5300, 7700, 2500),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH850": {"SO58": (5800, 10350, 1500),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH900": {"SO58": (7500, 10450, 3000),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH680": {"SY47": (6450, 7350, 7500),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH800": {"SY47": (7500, 8600, 7000),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH950": {"O58": (8850, 10000, 5500),
                         "DEFAULT": (3000, 11000, 2500)}
              }
_GRISM_NAMES = (
    ("SCFCGREL01", "SCFCGRLD01", "SCFCGRMB01", "SCFCGRMR01", "SCFCGRHDEC", "SCFCGRHD45", "SCFCGRHD52", "SCFCGRHD65",
     "SCFCGRHD68", "SCFCGRHD80", "SCFCGRHD95"),
    ("75", "150", "300B", "300R", "Echelle", "VPH450", "VPH520", "VPH650", "VPH680", "VPH800", "VPH950")
)
GRISM_NAMES = {}

for idx in range(len(_GRISM_NAMES[0])):
    GRISM_NAMES[_GRISM_NAMES[0][idx]] = _GRISM_NAMES[1][idx]

FILTER_MAP = dict(SCFCFLBU01='U',
                  SCFCFLBB01='B',
                  SCFCFLBV01='V',
                  SCFCFLBR01='R',
                  SCFCFLBI01='I',
                  SCFCFLN373='N373',
                  SCFCFLN386='N386',
                  SCFCFLN487='N487',
                  SCFCFLN502='N502',
                  SCFCFLN512='N512',
                  SCFCFLN642='N642',
                  SCFCFLN670='N670',
                  SCFCFLBSZ1='SDSS_z',
                  SCFCFLSO58='O58',
                  SCFCFLSY47='Y47',
                  SCFCFLL600='L600',
                  SCFCFLL550='L550',
                  SCFCFLLC50='C50'
                  )


def augment_filter_lookup_table(bandpass_database):
    """
    Add some filters to the internal database, based on information from Subaru website.
    """

    energy_bouds = dict(O58=(580 * units.nm, 1000 * units.nm),
                        Y47=(470 * units.nm, 910 * units.nm),
                        SDSS_z=(813.850 * units.nm, 1026.852 * units.nm),
                        L600=(370 * units.nm, 600 * units.nm),
                        L550=(340 * units.nm, 525 * units.nm),
                        C50=(500 * units.nm, 11000 * units.nm),
                        Hcont=((1.573 - 0.020 / 2) * units.um, (1.573 + 0.020 / 2.0) * units.um),
                        Brgamma=((2.166 - 0.032 / 2.) * units.um, (2.166 + 0.032 / 2) * units.um))

    for key in energy_bouds:
        bandpass_database.add_static_filter(key, energy_bouds[key])
