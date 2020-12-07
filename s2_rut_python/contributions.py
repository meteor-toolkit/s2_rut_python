# From:
# Gorrono et al (2018). Providing uncertainty estimates of the Sentinel-2 top-of-atmosphere measurements for
# radiometric validation activities. European Journal of Remote Sensing, 51(1), 650-666.
# https://doi.org/10.1080/22797254.2018.1471739

U_INDEPENDENT_CONTRIBUTIONS = ["Instrument_noise",
                               "ADC_quantisation",
                               "L1C_image_quantisation"]

U_COMMON_CONTRIBUTIONS = ["Crosstalk",
                          "Diffuser-absolute_knowledge",
                          "Diffuser-cosine_effect",
                          "Diffuser-straylight_residual",
                          "Diffuser-temporal_knowledge",
                          "DS_stability",
                          "Gamma_knowledge",
                          "OOF_straylight-random",
                          "OOF_straylight-systematic"]
