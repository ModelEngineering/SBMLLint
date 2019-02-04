"""Checks for static errors in a model."""

MODEL_TYPE_FILE = "model_type_file"
MODEL_TYPE_ANTIMONY = "model_type_antimony"
MODEL_TYPE_URL = "model_type_url"
MODEL_TYPE_BIOMODELS = "model_type_biomodels"


def checker(model,
    model_type=MODEL_TYPE_ANTIMONY,
    mass_balance_check="structured_names"):
  """
  Reports on errors found in a model
  :param object model: reaction network model
  :param str model_type: type of model used
  :param str mass_balance_check: how check for mass balance
  """
  pass
