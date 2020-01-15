import configparser
import os


def setup():
    # use config parser to get settings from config.ini
    settings_file = os.path.join(os.getcwd(), 'config.ini')
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read(settings_file)
    return settings
