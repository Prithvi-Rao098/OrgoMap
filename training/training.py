'''
API:
- Diagrams: https://coconut.naturalproducts.net/search
- List of datasets: https://github.com/kjappelbaum/awesome-chemistry-datasets
- Reactions: https://github.com/MolecularAI/Chemformer, https://github.com/doyle-lab-ucla/ochem-data
- Common compounds: https://github.com/cheminfo/pubchem
- "diagrams" API: https://huggingface.co/datasets/camel-ai/chemistry
'''

# pip install google-cloud-aiplatform
# gs://your-bucket/path/to/training_data.jsonl  
# gs://your-bucket/path/to/validation_data.jsonl  


from vertexai.preview.generative_models import SupervisedTuningJob
from vertexai import init

# Initialize Vertex AI
init(project="your-gcp-project-id", location="us-central1")

# Define file locations in your GCS bucket
TRAIN_URI = "gs://your-bucket/path/to/training_data.jsonl"
VAL_URI = "gs://your-bucket/path/to/validation_data.jsonl"

# Create and start the fine-tuning job
job = SupervisedTuningJob.create(
    display_name="gemini-flash-organicchem-tuning",
    model="gemini-1.5-flash",  # or gemini-1.0-flash
    training_data_uri=TRAIN_URI,
    validation_data_uri=VAL_URI,
    tuned_model_display_name="gemini-flash-chemistry-expert"
)

from vertexai.preview.generative_models import GenerativeModel

model = GenerativeModel(model_name="tunedModels/gemini-flash-chemistry-expert")

response = model.generate_content("What is the major product of 2-butene + HBr?")
print(response.text)
