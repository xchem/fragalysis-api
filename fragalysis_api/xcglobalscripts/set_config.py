import configparser
import os


def ConfigSetup():
    # use config parser to get settings from config.ini
    settings_file = os.path.join(os.path.dirname(__file__), "config.ini")
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read(settings_file)
    return settings
